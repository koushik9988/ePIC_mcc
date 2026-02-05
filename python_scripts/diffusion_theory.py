"""
Simple python script to analyze Electric field to study turbulence.
@kaushik kalita
@Date : 16/10/2025
Works with datasets of 1st paper/work.
"""
import numpy as np
import h5py
import sys
from os.path import join as pjoin
from scipy.constants import value as constants 
from scipy.stats import linregress  
import os
from os.path import join as pjoin
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


#----------------constanst--------------------

eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#----------------------------------------------

if len(sys.argv) != 3:
    print("Usage: python Dth_auto_full.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
file_name = 'result.h5'

output_dir = 'turbulence_plots'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

#open hdf file
f = h5py.File(pjoin(path, file_name), 'r')

#---------------------Read Metadata-----------------------------------------
metadata_group = f['/metadata']

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

species_metadata_path = f'/metadata_species/{particle_type}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']

#----------------------------------------------

#-------------------------norm------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)


L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi

DT = DT_coeff * (1.0 / W)
efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12


# ------------------------------------------------------
# Load full Electric Field data (time × space)
# ------------------------------------------------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}'][:]
    electric_field_data.append(EF_data)

EF_full = np.vstack(electric_field_data)   # shape: (Nt_full, Nx)
Nt_full, Nx = EF_full.shape

x  = np.linspace(0, NC, Nx, endpoint=False)
dx = x[1] - x[0]

dt_full = DT_coeff * write_interval * mfactor   # physical time step
time_full = np.arange(Nt_full) * dt_full


id1 = 1500   # in normalized time units; choose window after saturation
id2 = 4800

dt_field = dt_full   # same as above

start_idx = int(round(id1 / dt_field))
end_idx   = int(round(id2 / dt_field))

if end_idx > Nt_full:
    end_idx = Nt_full

print("Time window for Dth: t_norm in [", id1, ",", id2, "]")
print("start_idx:", start_idx, "end_idx:", end_idx)

EF = EF_full[start_idx:end_idx, :]
Nt = EF.shape[0]

kfft   = 2.0 * np.pi * np.fft.fftfreq(Nx, dx)
k_modes = kfft[:Nx // 2]                     # positive k
F_k     = np.fft.fft(EF, axis=1, norm='ortho')[:, :Nx // 2]

# power spectrum in k (time averaged, sliced window)
power_spectrum = np.mean(np.abs(F_k)**2, axis=0)

# detect dominant spectral peaks in k-space
peaks, _ = find_peaks(power_spectrum,height=np.max(power_spectrum) * 0.01,distance=2)

if len(peaks) == 0:
    print("No spectral peaks found in power_spectrum. Exiting.")
    sys.exit(0)

max_idx = np.argmax(power_spectrum)
k_val_norm = abs(k_modes[max_idx])
k_val      = k_val_norm / LDi   # physical k

print(f"Detected most unstable mode index: {max_idx}")
print("k_max (normalized k*LDi):", k_val * LDi)

# ---------------- Single-k Dth for reference (using k_max) ----------------
Ek_single = np.abs(F_k[:, max_idx]) * efield_norm * (2.0 / np.sqrt(NC))
Ek_max_single  = np.max(Ek_single)
Ek2_max_single = np.max(Ek_single**2)

print("Max field amplitude at k_max [V/m]:", Ek_max_single)

delta_vk_half_single = 2.0 * np.sqrt((e * Ek2_max_single) / (species_mass * k_val))
delta_vk_single      = 2.0 * delta_vk_half_single

print("Resonance velocity width (single k_max, normalized by LDi*wpi):",
      delta_vk_single / (LDi * wpi))

psdk_single = Ek_max_single * Ek_max_single

Dth_single = (np.pi * e**2 * psdk_single) / (2.0 * species_mass**2 * k_val * delta_vk_single)
Dth_single /= (wpi * (LDi * wpi)**2)

print("Theoretical Diffusion Coefficient Dth (single k_max):", Dth_single)


print("\nDetected dominant k-indexes (peaks):", peaks)

valid_indices = []
for idx in peaks:
    if abs(k_modes[idx]) != 0.0:
        valid_indices.append(idx)

Dth_all    = []
k_phys_all = []

for idx in valid_indices:

    k_norm = abs(k_modes[idx])
    k_phys = k_norm / L
    k_phys_all.append(k_phys)

    # E_k(t) amplitude in physical units
    Ek_t_amp = np.abs(F_k[:, idx]) * efield_norm * (2.0 / np.sqrt(NC))
    Ek_max   = np.max(Ek_t_amp)
    Ek2_max  = Ek_max**2

    ###
    #Ek_rms = np.sqrt(np.mean(Ek_t_amp**2))
    #Ek2_rms = Ek_rms**2
    #Ek2_max = Ek2_rms

    ###

    # resonance velocity width Δv_k (physical)
    delta_vk_half = 2.0 * np.sqrt((e * Ek2_max) / (species_mass * k_phys))
    delta_vk      = 2.0 * delta_vk_half

    # theoretical D(k) – same formula as single-k
    psdk = Ek2_max
    D_k = (np.pi * e**2 * psdk) / (2.0 * species_mass**2 * k_phys * delta_vk)
    D_k /= (W * (L * W)**2)

    Dth_all.append(D_k)

Dth_all    = np.array(Dth_all)
k_phys_all = np.array(k_phys_all)

# ------------------------------------------------------
# PHASE VELOCITY: from full EF history (not sliced)
# ------------------------------------------------------
# Spatial FFT on full EF
F_k_full = np.fft.fft(EF_full, axis=1, norm='ortho')[:, :Nx // 2]

vphi_full = []   # normalized v_phi = ω_norm / k_norm = v_phys / (LDi*wpi)

for idx in valid_indices:

    k_dim = k_modes[idx]    # normalized k
    if k_dim == 0.0:
        continue

    # mode time series over full simulation
    Ek_full = F_k_full[:, idx]

    # time FFT → ω(k)
    wspec_full = np.fft.fft(Ek_full)
    freqs_full = 2.0 * np.pi * np.fft.fftfreq(Nt_full, dt_full)   # ω in normalized units

    pos_mask = freqs_full > 0.0
    if not np.any(pos_mask):
        # fallback: use half spectrum if something odd
        pos_mask = slice(0, Nt_full // 2)

    abs_wspec = np.abs(wspec_full[pos_mask])
    omega_peak = freqs_full[pos_mask][np.argmax(abs_wspec)]

    # normalized phase velocity v_phi_norm = ω_norm / k_norm
    vphi_k = omega_peak / k_dim
    vphi_full.append(vphi_k)

vphi_full = np.array(vphi_full)

print("\n======= Phase Velocity extracted from COMPLETE EF history =======")
for i, (k_norm_val, vph_val, Dval) in enumerate(zip(k_phys_all * LDi, vphi_full, Dth_all)):
    print(f"Mode {i}: k*LDi = {k_norm_val:7.3f}   vphi_norm = {vph_val:9.3f}   D_th = {Dval:.4e}")

# ------------------------------------------------------
# VISUALIZATION: D_th(k) vs v_phi (vertical lines)
# ------------------------------------------------------
plt.figure(figsize=(8, 5))

for v, D in zip(vphi_full, Dth_all):
    plt.vlines(v, 0.0, D, lw=2.0)
    plt.scatter(v, D, s=35, color='red')

plt.xlabel(r"$v_\phi = \omega_k / k$")
#plt.xlim([5,20])
plt.ylabel(r"$D_{\mathrm{th}}$")
plt.title("Dth vs v_phase")
plt.grid(alpha=0.3)
plt.tight_layout()

fig_name = pjoin(path_fig, "Dth_vs_vphase.png")
plt.savefig(fig_name, dpi=300)
plt.show()
