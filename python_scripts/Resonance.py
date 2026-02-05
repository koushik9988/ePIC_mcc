import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks

# ------------------------------
# Usage Check
# ------------------------------
if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data_directory>")
    sys.exit(1)

path = sys.argv[1]
file_name = 'result.h5'

# Load HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')

# Metadata
metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
normscheme = metadata.attrs['norm_scheme']
density = metadata.attrs['density']


mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)


DT = DT_coeff * (1.0 / wpi)
dt_field = DT_coeff * write_interval * mfactor


efield_norm = (density * 1.602176565e-19 * LDi) / 8.85418782E-12 #for ion scale 

# Load Electric Field Data

electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}'] #* efield_norm
    electric_field_data.append(EF_data[:])

f.close()

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
time_full = np.arange(Nt) * dt_field

# Find dominat K mode via Spatial FFT
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)
k_modes = kfft[:Nx // 2]
F_k_pos = F_k[:, :Nx // 2]

power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)

print("Detected harmonic peak indices:", peaks)
print("Corresponding k values:", k_modes[peaks])

max_idx = peaks[np.argmax(power_spectrum[peaks])]
k_dominant = k_modes[max_idx]
print(f"Dominant k = {k_dominant:.4f}")

# Temporal FFT for Dominant k
E_k_t = F_k_pos[:, max_idx]
dt_snapshot = time_full[1] - time_full[0]

omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)
omega_pos_freqs = omega_freqs[:Nt // 2]
E_omega2 = np.abs(omega_spectrum[:Nt // 2])**2


# Find FWHM (Resonance Broadening)
peak_idx = np.argmax(E_omega2)
peak_val = E_omega2[peak_idx]
half_max = 0.5 * peak_val

# Find ω1 and ω2 around the main peak
left_idx = np.where(E_omega2[:peak_idx] <= half_max)[0]
right_idx = np.where(E_omega2[peak_idx:] <= half_max)[0]

if len(left_idx) > 0:
    omega1 = omega_pos_freqs[left_idx[-1]]
else:
    omega1 = omega_pos_freqs[0]

if len(right_idx) > 0:
    omega2 = omega_pos_freqs[peak_idx + right_idx[0]]
else:
    omega2 = omega_pos_freqs[-1]

fwhm = omega2 - omega1

print(f"Dominant ω = {omega_pos_freqs[peak_idx]:.4f}")
print(f"Resonance broadening (Δω = ω2 - ω1) = {fwhm:.4f}")

# Plot E(ω)^2 vs ω
plt.figure()
plt.plot(omega_pos_freqs, E_omega2, 'b-', lw=1.5)
plt.axvline(omega_pos_freqs[peak_idx], color='r', ls='--', label='Peak ω')
plt.axvline(omega1, color='g', ls='--', label='ω₁ (FWHM)')
plt.axvline(omega2, color='g', ls='--', label='ω₂ (FWHM)')
plt.title(f"E(ω)² vs ω for k = {k_dominant:.4f}\nΔω = {fwhm:.4f}")
plt.xlabel("ω")
plt.ylabel("|E(ω)|²")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
