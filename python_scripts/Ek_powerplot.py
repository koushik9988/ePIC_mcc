"""
This script plots the amplitude of the dominant wave mode vs time
and fits its exponential growth or damping rate.
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin


if len(sys.argv) != 4:
    print("Usage: python script.py <path_to_data> <fit_start_fraction> <fit_end_fraction>")
    sys.exit(1)

path = sys.argv[1]

print("Invalid fit range. fit1 and fit2 should be fractions (0 < fit1 < fit2 < 1).")
fit1 = float(sys.argv[2])   # Fraction of total time for start of fit
fit2 = float(sys.argv[3])   # Fraction of total time for end of fit

file_name = 'result.h5'

# ----------------------------- Read Metadata -----------------------------
f = h5py.File(pjoin(path, file_name), 'r')

metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']

iontemp = metadata_ion.attrs['temperature']
spwt_ion = metadata_ion.attrs['spwt']
ion_num = metadata_ion.attrs['num_particles']
ion_mass = metadata_ion.attrs['mass']
ion_charge = metadata_ion.attrs['charge']

metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
write_interval_phase = metadata.attrs['write_int_phase']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
normscheme = metadata.attrs['norm_scheme']
density = metadata.attrs['density']


mfactor = wpi / wpe
if normscheme in [1, 2]:
    mfactor = 1

L = LDe
if normscheme in [2,4]:
    L = LDi

print("mfactor:", mfactor)
print("normalization scheme:", normscheme)

efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12

# ----------------------------- Load Electric Field Data -----------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}'][:]
    electric_field_data.append(EF_data)

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
x = np.linspace(0, NC, EF.shape[1], endpoint=False)
dx = x[1] - x[0]

# ----------------------------- Spatial FFT -----------------------------
F_k = np.fft.fft(EF, axis=1, norm='ortho') * efield_norm * (2 / np.sqrt(NC))
print("fk multi factor ", efield_norm * (2 / np.sqrt(NC)))
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)

# ----------------------------- Time and Fit Window ---------------------
Nt = EF.shape[0]
time_full = np.arange(Nt) * DT_coeff * write_interval * mfactor

# Ensure correct integer indices for slicing
t1 = int(fit1 * Nt)
t2 = int(fit2 * Nt)

time_fit = time_full[t1:t2]

power_spectrum = np.mean(np.abs(F_k)**2, axis=0)
k_max_idx = np.argmax(power_spectrum[:NC // 2])
k_max = kfft[k_max_idx]

E_k_max_full = np.abs(F_k[:, k_max_idx])
E_k_max_fit = E_k_max_full[t1:t2]

# Filter valid values
valid = E_k_max_fit > 0
time_valid = time_fit[valid]
log_E = np.log(E_k_max_fit[valid])

# Fit exponential (ln|E| = γt + c)
coeffs = np.polyfit(time_valid, log_E, 1)
gamma = coeffs[0]
gamma_norm = gamma #/ wpi

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


ax2.semilogy(time_full, E_k_max_full, label='PIC Data', color='blue')
ax2.semilogy(time_fit, np.exp(coeffs[1] + coeffs[0] * time_fit),'r--', label=f'Fit: $\\gamma/\\omega_{{pe}}$ = {gamma_norm:.4f}')
ax2.set_xlabel(r'$\omega_{pi}t$')
ax2.set_ylabel(r'$|E(k_{max})|$')
ax2.set_title('Growth or Damping of Dominant Mode')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig(pjoin(path, 'growth_analysis.pdf'), dpi=300)
plt.show()

print(f"Dominant wavenumber (k_max): {k_max:.4f}")
print(r"Growth rate $\gamma$ (normalized by ω_pe): {gamma_norm:.4f}")

#k_launched = (5 * 2 * np.pi) / NC
#k_lambdaD = k_launched * LDe
#print(f"k value of launched wave: {k_launched:.4f}")
#print(f"kλ_D (launched): {k_lambdaD:.4f}")
#print(f"k_max/k_launched = {k_max/k_launched:.4f}")
#rint(f"k_launched/k_max = {k_launched/k_max:.4f}")
