import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks


if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data_directory>")
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
normscheme = metadata.attrs['norm_scheme']


mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)


DT = DT_coeff * (1.0 / wpi)

# Load Electric Field Data
electric_field_data = []

time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

f.close()

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx, endpoint=False)
dx = x[1] - x[0]
time_full = np.arange(Nt) * DT_coeff * write_interval*mfactor #normalized time

# Spatial FFT and Peak Finding
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx) #normalized k we havenot used dx = debye lenght so k is normalized
k_modes = kfft[:Nx // 2] # Positive k modes
F_k_pos = F_k[:, :Nx // 2] # Positive k modes

power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)

    
print("Detected harmonic peak indices:", peaks)
print("Corresponding k values:", k_modes[peaks])
max_idx = peaks[np.argmax(power_spectrum[peaks])]

#Isolate the dominant mode's wavenumber and time-series data
k_dominant = k_modes[max_idx]
E_k_t = F_k_pos[:, max_idx]  #Complex amplitude E(k, t)


dt_snapshot = time_full[1] - time_full[0]
omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)

omega_pos_freqs = omega_freqs[:Nt // 2] # Positive frequencies
omega_pos_spectrum = omega_spectrum[:Nt // 2] # Positive frequencies modes
omega_peak_idx = np.argmax(np.abs(omega_pos_spectrum))
omega_r = omega_pos_freqs[omega_peak_idx]

v_phase = omega_r / k_dominant #normalized phase velocity

print(f"Dominant Wavenumber (k): {k_dominant:.4f}")
print(f"Dominant Frequency (ω_r): {omega_r:.4f}")
print(f"Phase Velocity (v_φ = ω_r/k * vthi): {v_phase:.4f}")
