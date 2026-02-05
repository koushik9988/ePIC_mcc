import numpy as np 
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks


if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]


file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

# Read Metadata 
metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
density = metadata.attrs['density']
normscheme = metadata.attrs['norm_scheme']

DT = DT_coeff * (1.0 / wpe)

mFactor = wpi/wpe


mFactor = wpi/wpe 
if normscheme in [1, 2]:
    mFactor = 1



L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi


efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12

#  Load Electric Field Data 
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    #EF_data = f[f'fielddata/den_electron/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
x = np.linspace(0, NC, EF.shape[1], endpoint=False)
dx = x[1] - x[0]


time_full = np.arange(EF.shape[0]) * DT_coeff * write_interval*mFactor

# Spatial FFT 
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)
k_modes = kfft[:NC // 2]
F_k_pos = F_k[:, :NC // 2] * efield_norm * (2/np.sqrt(NC))  # cut half spectrum (Positive k modes only)

# Peak finding
# distance = minimum number of samples between neighboring peaks.
#height = minimum height threshold for peaks to be detected.
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.00, distance = 2)

print("Detected harmonic peak indices:", peaks)
print("Corresponding k values:", k_modes[peaks])

#-- harmonics numbers
k_peaks = k_modes[peaks]
k_fundamental = k_peaks[0] # First harmonic

# Calculate the ratio
harmonic_ratios = k_peaks / k_fundamental

print(" (k_peak / k_fundamental):", harmonic_ratios)
print("Approximate harmonic numbers:", np.round(harmonic_ratios).astype(int))

#----------------------------

# ----------------------------- Plot Power Spectrum with Peak Markers -----------------------------
fig_ps, ax_ps = plt.subplots()

ax_ps.plot(k_modes, power_spectrum, label='Power Spectrum')
#ax_ps.semilogy()
ax_ps.plot(k_modes[peaks], power_spectrum[peaks], 'ro', label='detected modes')

# maximum peak
max_idx = peaks[np.argmax(power_spectrum[peaks])]
ax_ps.plot(k_modes[max_idx], power_spectrum[max_idx], 'ko', markersize=8,label=f'Max Peak at k ≈ {k_modes[max_idx]:.2f}')
ax_ps.set_xlabel('$k$')
ax_ps.set_ylabel('$|E(k)|^2$')
ax_ps.grid(True)
ax_ps.legend()
plt.tight_layout()
plt.savefig(pjoin(path, 'power_spectrum_peaks.pdf'), dpi=300)

# Energy Evolution of Detected Harmonics -----------------------------
energy_spectrum = np.abs(F_k_pos)**2 

fig, ax = plt.subplots(figsize=(10, 6))
#for idx in peaks:
#    energy = energy_spectrum[:, idx]
#    ax.plot(time_full, energy, label=f'$k \\approx {k_modes[idx]:.2f}$')

num_peaks_to_plot = min(5, len(peaks))
for idx in peaks[:num_peaks_to_plot]:
#for idx in peaks:
    energy = energy_spectrum[:, idx]
    ax.plot(time_full, np.log(energy), label=f'$k \\approx {k_modes[idx]:.2f}$')


ax.set_xlabel('Time')
ax.set_ylabel('$|E_k|^2$')
#ax.set_yscale('log')
ax.set_title('Time Evolution of Harmonics')
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(pjoin(path, 'energy_exchange_between_modes.pdf'), dpi=300)
plt.show()
