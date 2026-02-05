import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os

from scipy.interpolate import interp1d


if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data_directory>")
    sys.exit(1)

path = sys.argv[1]
file_name = 'result.h5'

plot_path = './turbulence'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

f = h5py.File(pjoin(path, file_name), 'r')

metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']


iontemp = metadata_ion.attrs['temperature']
spwt_ion = metadata_ion.attrs['spwt']
ion_num = metadata_ion.attrs['num_particles']
ion_mass = metadata_ion.attrs['mass']
ion_charge = metadata_ion.attrs['charge']

# Read Metadata
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

mfactor = wpi/wpe 
if normscheme in [1, 2]:
    mfactor = 1

L = LDe
if normscheme == 2:
    L = LDi

print("mfactor:", mfactor)
print("normalization scheme:", normscheme)

efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12
#load field
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for ts in time_steps:
    EF_data = f[f'fielddata/efield/{ts}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx, endpoint=False)
dx = x[1] - x[0]

print(f"Data shape: EF({Nt}, {Nx})")

# Spatial FFT
F_k = np.fft.fft(EF, axis=1, norm='ortho') * efield_norm * (2 / np.sqrt(NC))
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)

# Temporal FFT
F_k_t = np.fft.fft(F_k, axis=0, norm='ortho')
omega_fft = 2 * np.pi * np.fft.fftfreq(Nt, DT_coeff * write_interval * mfactor)

# Shift FFT outputs to center zero frequency/wavenumber
F_k_t_shifted = np.fft.fftshift(F_k_t, axes=(0, 1))
k_shifted = np.fft.fftshift(kfft)
omega_shifted = np.fft.fftshift(omega_fft)

# Compute power spectrum
dispersion_power = np.abs(F_k_t_shifted)**2
dispersion_power /= np.max(dispersion_power)  # normalize


power_log = np.log10(dispersion_power)

power_log = power_log[:, ::-1]
k_shifted = -k_shifted[::-1]



################## fidn domiant w-k pair ######################################
psd_vs_k = np.mean(dispersion_power, axis=0)  # Shape (Nx,)

k_peak_indices, _ = find_peaks(psd_vs_k, height=np.max(psd_vs_k) * 0.01, distance=2)
dominant_k_values = -k_shifted[k_peak_indices]

dominant_omega_values = []
for k_idx in k_peak_indices:
    psd_vs_omega_at_k = dispersion_power[:, k_idx]  # Shape (Nt,)
    
    omega_peak_indices, _ = find_peaks(psd_vs_omega_at_k, height=np.max(psd_vs_omega_at_k) * 0.01, distance=1)
    
    if len(omega_peak_indices) > 0:
        strongest_omega_idx_in_peaks = np.argmax(psd_vs_omega_at_k[omega_peak_indices])
        strongest_omega_idx = omega_peak_indices[strongest_omega_idx_in_peaks]
        dominant_omega_values.append(omega_shifted[strongest_omega_idx])
    else:
        dominant_omega_values.append(np.nan)

dominant_omega_values = np.array(dominant_omega_values)

valid_indices = ~np.isnan(dominant_omega_values)
dominant_k_plot = dominant_k_values[valid_indices]
dominant_omega_plot = dominant_omega_values[valid_indices]
print(f"Found {len(dominant_k_plot)} valid (k, ω) pairs to plot.")

for i in range(len(dominant_k_plot)):
    print(dominant_k_plot[i],":",dominant_omega_plot[i])
#################


# full (ω–k) spectrum
fig, ax = plt.subplots(figsize=(7, 5))
extent = [k_shifted[0], k_shifted[-1], omega_shifted[0], omega_shifted[-1]]

im = ax.imshow(power_log,extent=extent, origin='lower', aspect='auto', interpolation="bilinear", cmap='plasma', vmin=-10, vmax=0)

#ax.scatter(dominant_k_plot, dominant_omega_plot, color='cyan', s=25, marker='x', label='Dominant Modes')
#ax.scatter(dominant_k, dominant_omega, color='cyan', s=50, label='Dominant (k,ω)')

ax.set_xlabel(r"$k$")
ax.set_ylabel(r"$\omega$")
#ax.set_xlim(-2,2)
#ax.set_ylim(-10,10)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$\log_{10}|E(k,\omega)|^2$")
plt.tight_layout()
print(pjoin(path_fig, "dispersion_psd.png"))
plt.savefig(pjoin(path_fig, "dispersion_psd1.png"), dpi=300)
plt.show()
