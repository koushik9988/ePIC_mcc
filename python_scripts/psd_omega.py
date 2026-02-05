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

# ----------------------------- Read Metadata -----------------------------
metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
LDe = metadata.attrs['LDe']

DT = DT_coeff * (1.0 / wpe)

# ----------------------------- Load Electric Field Data -----------------------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
x = np.linspace(0, NC, EF.shape[1], endpoint=False)
dx = x[1] - x[0]
Nt = EF.shape[0]
time_full = np.arange(EF.shape[0]) * DT_coeff * write_interval
dt_eff = DT_coeff * write_interval  # Time step in units of 1/ω_pe

# ----------------------------- Temporal Power Spectrum |E(omega)|^2 -----------------------------
F_omega = np.fft.fft(EF, axis=0, norm='ortho')  # FFT along time axis
power_omega = 2 * np.abs(F_omega)**2 / (Nt * dt_eff)  # PSD with proper normalization
total_power_vs_omega = np.sum(power_omega, axis=1)  # Sum over spatial points

# Frequency axis (angular frequency in ω/ω_pe)
freqs = 2 * np.pi * np.fft.fftfreq(Nt, d=dt_eff)  # Angular frequencies (already in ω/ω_pe)
omega = np.fft.fftshift(freqs)  # Shift for plotting
omega_normalized = omega  # No further normalization needed
total_power_vs_omega_shifted = np.fft.fftshift(total_power_vs_omega)

# ----------------------------- Peak Detection -----------------------------
# Detect peaks in the positive frequency range
positive_omega_indices = np.where(omega_normalized >= 0)[0]
positive_omega = omega_normalized[positive_omega_indices]
positive_power = total_power_vs_omega_shifted[positive_omega_indices]
peaks, _ = find_peaks(positive_power, height=np.max(positive_power) * 0.04, distance=2)

print("Detected frequency peak indices (in positive omega):", peaks)
print("Corresponding ω/ω_pe:", positive_omega[peaks])
if len(peaks) > 0:
    f_peaks = positive_omega[peaks]
    f_fundamental = f_peaks[0]  # First peak as fundamental
    frequency_ratios = f_peaks / f_fundamental
    print("Frequency ratios (ω_peak / ω_fundamental):", frequency_ratios)
    print("Approximate harmonic numbers:", np.round(frequency_ratios).astype(int))
else:
    print("No peaks detected in the frequency spectrum.")

# ----------------------------- Plot |E(omega)|^2 vs omega -----------------------------
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(omega_normalized, total_power_vs_omega_shifted, color='darkblue', label='Power Spectrum')
if len(peaks) > 0:
    ax.plot(positive_omega[peaks], positive_power[peaks], 'ro', label='Detected modes')
    max_idx = peaks[np.argmax(positive_power[peaks])]
    ax.plot(positive_omega[max_idx], positive_power[max_idx], 'ko', markersize=8,
            label=f'Max Peak at ω/ω_pe ≈ {positive_omega[max_idx]:.2f}')
ax.set_xlabel(r'$\omega / \omega_{pe}$')
ax.set_ylabel(r'Power Spectral Density $|E(\omega)|^2$')
ax.set_title('Field Energy Spectrum vs Frequency')
ax.set_xlim(0, max(omega_normalized))  # Plot only positive omega
ax.set_yscale('log')  # Log scale for wide dynamic range
ax.grid(True, which='both', linestyle='--', alpha=0.6)
ax.legend()
plt.tight_layout()
plt.savefig(pjoin(path, 'field_energy_vs_omega.pdf'), dpi=300)
plt.show()

f.close()