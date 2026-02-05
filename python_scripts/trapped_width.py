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

e = np.abs(ion_charge)

mfactor = wpi/wpe 
if normscheme in [1, 2]:
    mfactor = 1

print("mfactor:", mfactor)
print("normalization scheme:", normscheme)

efield_norm = (density * 1.602176565e-19 * LDi) / 8.85418782E-12

# Load Electric Field Data
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])
EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
time_full = np.arange(Nt) * DT_coeff * write_interval * mfactor  # normalized time
dt_field = DT_coeff * write_interval * mfactor
dt_snapshot = dt_field #

# Spatial FFT and Peak Finding
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)  # normalized k (dimensionless)
k_modes = kfft[:Nx // 2]  # Positive k modes
F_k_pos = F_k[:, :Nx // 2]  # Positive k modes
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)
print("Detected harmonic peak indices:", peaks)
print("Corresponding k values:", k_modes[peaks])

max_idx = peaks[np.argmax(power_spectrum[peaks])]
k_dominant = k_modes[max_idx]
E_k_t = F_k_pos[:, max_idx]  # Amplitude E(k, t) (complex)
k_phys = k_dominant / LDi  # physical wavenumber (1/m)

# Temporal FFT for frequency
omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)
omega_pos_freqs = omega_freqs[:Nt // 2]  # Positive frequencies
omega_pos_spectrum = omega_spectrum[:Nt // 2]  # Positive frequencies
omega_peak_idx = np.argmax(np.abs(omega_pos_spectrum))
omega_r = omega_pos_freqs[omega_peak_idx]
v_phase = omega_r / k_dominant  # phase velocity in m/s

print(f"Dominant Wavenumber (k_dimless): {k_dominant:.6e}")
print(f"Dominant Wavenumber (k_phys, 1/m): {k_phys:.6e}")
print(f"Dominant Frequency (ω_r, rad/s): {omega_r:.6e}")
print(f"Phase Velocity (v_φ = ω_r/k): {v_phase:.6e} m/s")


#.........................Bounce frequency calculation................................................
E_k_amp = np.abs(E_k_t) * efield_norm * (2/np.sqrt(NC)) # as norm = ortho does divide by 1/sqrt(N) so for complete scaling we need 1/N
"""
fft result are proportional to have to scale by 1/N or in case of ortho(where persaval idenity is preserved) 1/sqrt(N) 
and a multiplication factor of 2 comes because energy is divided into both positive and negative spectra , so positive spectra carry half and negative half.

"""
omega_b_avg = (np.sqrt(k_phys * np.mean(E_k_amp) * (e / ion_mass)))/wpi
omega_b_max = (np.sqrt(k_phys * np.max(E_k_amp) * (e/ ion_mass)))/wpi
print(f"Average Bounce Frequency (ω_b, rad/s): {omega_b_avg:.6e}")
print(f"Max Bounce Frequency (ω_b, rad/s): {omega_b_max:.6e}")

# Trapping width calculation
#E_k_mean = np.mean(E_k_amp)
E_k_max = np.max(E_k_amp)
#average after psd reach max value to end
max_idx_time = np.argmax(E_k_amp)
E_k_mean = np.mean(E_k_amp[max_idx_time:])



delta_v_mean = 2.0 * np.sqrt((e * E_k_mean) / (ion_mass * np.abs(k_phys))) #half width
delta_v_max =  2.0 * np.sqrt(2 * e * E_k_max / (ion_mass * np.abs(k_phys))) #half width
#delta_v_diff = np.sqrt(2 * e * (E2 - E1) / (mi * np.abs(k_phys))) #half width

delta_v_mean_norm = delta_v_mean / (wpi * LDi)
delta_v_max_norm = delta_v_max / (wpi * LDi)
#delta_v_diff_norm = delta_v_diff / (wpi * LDi)
 
print(f"Trapping width (Δv) using mean E_k: {delta_v_mean:.6e} m/s (normalized: {delta_v_mean_norm:.6e})")
print(f"Trapping width (Δv) using max E_k: {delta_v_max:.6e} m/s (normalized: {delta_v_max_norm:.6e})")
#print(f"Trapping width (Δv) using (E_k_max - E_k_min): {delta_v_diff:.6e} m/s (normalized: {delta_v_diff_norm:.6e})")

# Trapped velocity range in lab frame (normalized)
v_min_mean = (v_phase - delta_v_mean_norm)
v_max_mean = (v_phase + delta_v_mean_norm)
v_min_max = (v_phase - delta_v_max_norm)
v_max_max = (v_phase + delta_v_max_norm)


print(f"Trapped velocity range in lab frame (using mean E_k, normalized): from {v_min_mean:.6e} to {v_max_mean:.6e}")
print(f"Trapped velocity range in lab frame (using max E_k, normalized): from {v_min_max:.6e} to {v_max_max:.6e}")

