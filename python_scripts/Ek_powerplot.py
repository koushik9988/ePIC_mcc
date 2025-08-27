import numpy as np 
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin

if len(sys.argv) != 4:
    print("Usage: python script.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]

fit1 = float(sys.argv[2])
fit2 = float(sys.argv[3])

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

EF = np.vstack(electric_field_data)  # Shape: (time_steps, spatial_points)
x = np.linspace(0, NC, EF.shape[1], endpoint=False)
dx = x[1] - x[0]

# ----------------------------- Time Slicing -----------------------------
wpet_1 = NUM_TS * DT_coeff * fit1  # start time for fit
wpet_2 = NUM_TS * DT_coeff * fit2 # end time for fit

y1 = int(wpet_1 / (DT_coeff * write_interval))
y2 = int(wpet_2 / (DT_coeff * write_interval))

# ----------------------------- Spatial FFT -----------------------------
F_k = np.fft.fft(EF, axis=1, norm='ortho')  # Full time FFT
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)

# Time array for full and fitting interval
time_full = np.arange(EF.shape[0]) * DT * write_interval
time_fit = time_full[y1:y2]

# Get E(k_max)
power_spectrum = np.mean(np.abs(F_k)**2, axis=0)
k_max_idx = np.argmax(power_spectrum[:NC // 2])
k_max = kfft[k_max_idx]

E_k_max_full = np.abs(F_k[:, k_max_idx])
E_k_max_fit = E_k_max_full[y1:y2]

# Filter for valid (non-zero) values
valid = E_k_max_fit > 0
log_E = np.log(E_k_max_fit[valid])
time_valid = time_fit[valid]

# Fit exponential to log(|E|)
coeffs = np.polyfit(time_valid, log_E, 1)
gamma = coeffs[0]
gamma_norm = gamma / wpe

# ----------------------------- Plotting -----------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# --- Power Spectrum ---
ax1.plot(kfft[:NC // 2], power_spectrum[:NC // 2], color='blue')
ax1.axvline(k_max, color='red', linestyle='--', label=f'$k_{{max}}$ = {k_max:.3f}')
ax1.set_xlabel('$k$')
ax1.set_ylabel('$|E(k)|^2$')
ax1.set_title('Spatial Power Spectrum')
ax1.grid(True)
ax1.legend()

# --- Growth/Decay of Mode ---
ax2.semilogy(time_full, E_k_max_full, label='PIC Data', color='blue')
ax2.semilogy(time_fit, np.exp(coeffs[1] + coeffs[0]*time_fit), 'r--', label=f'Fit: $\\gamma/\\omega_{{pe}}$ = {gamma_norm:.4f}')
ax2.set_xlabel('Time [$\omega_{pe}^{-1}$]')
ax2.set_ylabel('$|E(k_{max})|$')
ax2.set_title('Growth or Damping of Dominant Mode')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig(pjoin(path, 'growth_analysis.pdf'), dpi=300)
plt.show()

# ----------------------------- --------------------------------
print(f"Dominant wavenumber (k_max): {k_max:.4f}")
print(f"Growth rate γ (normalized by ω_pe): {gamma_norm:.4f}")
k_launched = (5*2*np.pi)/NC
k_lambdaD = k_launched*LDe
print("k value of launched wave",k_launched)
print(k_max/k_launched,"ratio",k_launched/k_max)