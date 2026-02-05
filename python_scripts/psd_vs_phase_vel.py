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

mFactor = wpi/wpe 
if normscheme in [1, 2]:
    mFactor = 1

print("mfactor:", mFactor)
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
time_full = np.arange(Nt) * DT_coeff * write_interval * mFactor  # normalized time
dt_field = DT_coeff * write_interval * mFactor
dt_snapshot = dt_field #/ wpi  # seconds

############################################################################
# Spatial FFT → E(k, t)
E_k_t = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)
#k_pos = kfft > 0  # only positive k
#k_modes = kfft[k_pos]
#E_k_t = E_k_t[:, k_modes]

# Keep only positive k and skip k=0 mode
k_modes = kfft[1:Nx // 2]
E_k_t = E_k_t[:, 1:Nx // 2]

#k_modes = kfft[:Nx // 2]
#E_k_t = E_k_t[:, :Nx // 2]

# For each k, do temporal FFT to find ω(k)
Nt = EF.shape[0]        
omega_vals = 2 * np.pi * np.fft.fftfreq(Nt, d=dt_snapshot)
omega_pos_mask = omega_vals > 0
omega_vals = omega_vals[omega_pos_mask]

v_ph_list = []
E2_list = []

for ik, k in enumerate(k_modes):
    # Temporal FFT for this k
    E_omega = np.fft.fft(E_k_t[:, ik], norm='ortho')
    E_omega = E_omega[omega_pos_mask]
    P_omega = np.abs(E_omega)**2
    
    # Find dominant frequency for this k
    omega_k = omega_vals[np.argmax(P_omega)]

    v_ph = omega_k / k
   
    # Average E(k,t)^2 as overall field power for that mode
    E2 = np.mean(np.abs(E_k_t[:, ik])**2)* efield_norm**2
    
    v_ph_list.append(v_ph)
    E2_list.append(E2)

v_ph_array = np.array(v_ph_list)
E2_array = np.array(E2_list)

plt.figure()
#plt.scatter(v_ph_array, (E2_array / np.max(E2_array)), s=10, color='dodgerblue', edgecolors='k', alpha=0.8)
plt.scatter(v_ph_array, np.log(E2_array), s=10, color='dodgerblue', edgecolors='k', alpha=0.8)
plt.xlabel(r'$\omega_k / k$')
plt.ylabel('$|E(k)|^2$')
plt.xlim(0, 50)
plt.grid(False)
plt.tight_layout()
plt.savefig(pjoin(path_fig, 'E_k2_vs_vph.png'), dpi=200)
plt.show()

