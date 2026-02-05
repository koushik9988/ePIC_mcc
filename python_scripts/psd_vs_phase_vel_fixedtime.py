"""
@kaushik kalita
@Date : 29/10/2025
Script to analyze energy content in diffrent wave modes/ phase velcoities.
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation
from scipy.signal import find_peaks


# Argument
if len(sys.argv) != 2:
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
#particle_type = sys.argv[2]

plot_path = './turbulence'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
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

#species_metadata_path = f'/metadata_species/{particle_type}'
#metadata_species = f[species_metadata_path]

#species_mass = metadata_species.attrs['mass']
#species_charge = metadata_species.attrs['charge']
#species_temp = metadata_species.attrs['temperature']
#species_spwt = metadata_species.attrs['spwt']
#pecies_num = metadata_species.attrs['num_particles']

#----------------------------------------------

#-------------------------norm------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

DT = DT_coeff * (1.0 / wpi)

efield_norm = (density * 1.602176565e-19 * LDi) / 8.85418782E-12

# Load Electric Field Data
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:]) #shallow copy not reference copy
EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
time_full = np.arange(Nt) * DT_coeff * write_interval * mfactor  # normalized time
dt_field = DT_coeff * write_interval * mfactor
dt_snapshot = dt_field #/ wpi  # seconds

############################################################################
# Spatial FFT → E(k, t)
E_k_t = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, d=dx) # 2*pi*n/L ,L = Nx*dx and n= 0,1,2,...Nx-1
#k_pos = kfft > 0  # only positive k
#k_modes = kfft[k_pos]
#E_k_t = E_k_t[:, k_modes]

# Keep only positive k and skip k=0 mode
k_modes = kfft[1:Nx // 2]
E_k_t = E_k_t[:, 1:Nx // 2]

#k_modes = kfft[:Nx // 2]
#E_k_t = E_k_t[:, :Nx // 2]

# Compute dominant phase velocities from full time series
omega_vals = 2 * np.pi * np.fft.fftfreq(Nt, d=dt_snapshot)
omega_pos_mask = omega_vals > 0
omega_vals = omega_vals[omega_pos_mask]

v_ph_list = []

for ik, k in enumerate(k_modes): # enumerate() gives both the index (ik) and the value (k) at each loop.
    # Temporal FFT for this k
    E_omega = np.fft.fft(E_k_t[:, ik], norm='ortho')
    E_omega = E_omega[omega_pos_mask]
    P_omega = np.abs(E_omega)**2
    
    # Find dominant frequency for this k
    omega_k = omega_vals[np.argmax(P_omega)]

    v_ph = omega_k / k
   
    v_ph_list.append(v_ph)

v_ph_array = np.array(v_ph_list)


scaling = efield_norm**2 * (4 / NC)

# Plot interval 
plot_interval = 10
if Nt < plot_interval:
    plot_interval = 1

for t in range(0, Nt, plot_interval):
    # Instantaneous E2 at time t
    E2_t = np.abs(E_k_t[t, :])**2 * scaling
    E2_t_normalized = E2_t / np.max(E2_t) if np.max(E2_t) > 0 else E2_t

    plt.figure()
    plt.scatter(v_ph_array, E2_t_normalized, s=30, color='dodgerblue', edgecolors='k', alpha=0.8)
    plt.semilogy()
    #plt.scatter(v_ph_array, np.log(E2_array), s=30, color='dodgerblue', edgecolors='k', alpha=0.8)
    plt.xlabel(r'$\omega_k / k$')
    plt.ylabel('$|E(k,t)|^2$')
    plt.xlim(0, 50)
    plt.title(f'Time step {t}, t = {time_full[t]:.2f}')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(pjoin(path_fig, f'E_k2_vs_vph_t{time_full[t]:.2f}.png'), dpi=900)
    plt.close()