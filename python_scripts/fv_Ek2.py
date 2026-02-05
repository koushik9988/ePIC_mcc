"""
Simple python script to analyze Electric field to study turbulence.
@kaushik kalita
@Date : 11/11/2025
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
import matplotlib.animation as animation


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
species_name = sys.argv[2]
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

species_metadata_path = f'/metadata_species/{species_name}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']
species_spwt = metadata_species.attrs['spwt']

print(species_spwt)
#----------------------------------------------

#-------------------------norm------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

DT = DT_coeff * (1.0 / wpi)

efield_norm = (density * 1.602176565e-19 * LDi) / 8.85418782E-12


#--------load field data--------------------------------
electric_field_data = []
efield_time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for time_step in efield_time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])
EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape

x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
dt_field = DT_coeff * write_interval * mfactor
dt_snapshot = dt_field

# Spatial FFT → E(k, t)
E_k_t_all = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)

# Keep only positive k and skip k=0 mode
k_modes = kfft[1:Nx // 2]
E_k_t_all = E_k_t_all[:, 1:Nx // 2] # Shape (Nt, Nk)

# Compute dominant phase velocities from full time series
omega_vals = 2 * np.pi * np.fft.fftfreq(Nt, d=dt_snapshot)
omega_pos_mask = omega_vals > 0
omega_vals_pos = omega_vals[omega_pos_mask]

v_ph_list = []
for ik, k in enumerate(k_modes):
    # Temporal FFT for this k
    E_omega = np.fft.fft(E_k_t_all[:, ik], norm='ortho')
    E_omega = E_omega[omega_pos_mask]
    P_omega = np.abs(E_omega)**2
    
    # Find dominant frequency for this k
    omega_k = omega_vals_pos[np.argmax(P_omega)]
    v_ph = omega_k / k
    v_ph_list.append(v_ph)

v_ph_array = np.array(v_ph_list) # This is the static x-axis for the PSD plot
scaling = efield_norm**2 * (4 / NC)


#-----------Computing E-field power bounds...
E2_all_log = []
for t in range(0,len(efield_time_steps)):
    E2_t = np.abs(E_k_t_all[t, :])**2 * scaling
    E2_all_log.append(np.log(E2_t))

E2_all_log = np.array(E2_all_log)
E2_log_min = np.min(E2_all_log)
E2_log_max = np.max(E2_all_log)

print(f"Global log(|E|^2) range: {E2_log_min:.2f} to {E2_log_max:.2f}")
#-----------------------------------------------------------------------


#--------------------Find global velocity min/max-----------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1
# This is the time array corresponding to the particle snapshots
ts_particle = f["time_var/kinetic_energy"][:, 0] * mfactor


v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    v_min = min(v_min, np.min(data_species))
    v_max = max(v_max, np.max(data_species))
print(f"Velocity range found: {v_min:.2f} to {v_max:.2f}")

#-----------------Computing global f(v) max for normalization----------------
nbins = 100
fmax_global = 0.0

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    hist_species, _ = np.histogram(data_species, bins=nbins, range=(v_min, v_max), density=True)
    fmax_global = max(fmax_global, np.max(hist_species))

print(f"Global f(v) max value: {fmax_global:.3e}")

#-------------------------------------------------------------------------------

j0 = 0
data_species_0 = f[f"particle_{species_name}/{j0}"][:, 1]
nbins = 100
hist_species_0, bins = np.histogram(data_species_0, bins=nbins, range=(v_min, v_max), density=True)
bin_centers_0 = 0.5 * (bins[:-1] + bins[1:])


# Create figure with 2 subplots, sharing the x-axis
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)
plt.subplots_adjust(hspace=0.1) # Reduce space between plots

def animate(i):
    # i is the frame number (0, 1, 2, ...)
    # j is the simulation time step (0, 100, 200, ...)
    j = i * write_interval_phase
    
    ax1.clear()
    
    data_species = f[f"particle_{species_name}/{j}"][:, 1]

    hist_species, bins_hist = np.histogram(data_species, bins=nbins, range=(v_min, v_max), density=True)
    #hist_species *= species_spwt
    bin_centers = 0.5 * (bins_hist[:-1] + bins_hist[1:])
    if np.max(hist_species) > 0:
        hist_species /= fmax_global 

    #ax1.plot(bin_centers_0, hist_species_0, linestyle='--', color='gray', label='Initial')
    ax1.plot(bin_centers, hist_species, color='blue')
    ax1.set_ylim(0,1.05)
    ax1.set_ylabel('f(v)')
    #ax1.legend(loc='upper right')
    ax1.grid(True, linestyle=':', alpha=0.6)

    ax2.clear()
    
    #efield_time_steps = [0*write_interval,1*write_interval,2*write_interval,...........3001*write_interval] (size is 3001)
    #index(j) searches for value j in the efield_time_steps array and if it find then return the index of it (index is in betweeb 0 and 3001)
    t = efield_time_steps.index(j) 
    # E_k_t_all[t, :] is the E(k) spectrum at this time
    E2_t = np.abs(E_k_t_all[t, :])**2 * scaling
    E2_t_normalized = E2_t

    ax2.scatter(v_ph_array, np.log(E2_t_normalized), s=10, color='dodgerblue', edgecolors='k', alpha=0.8)
    ax2.set_ylabel('$|E(k)|^2$')
    #ax2.set_ylim(E2_log_min,(E2_log_max+0.5))
    ax2.set_xlabel('$\omega_k / k$')
    ax2.grid(True)

    # Set shared x-axis limits from particle data
    ax1.set_xlim(v_min, v_max)
    ax2.set_xlim(v_min, v_max)

    return ax1, ax2

def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS_PHASE - 1)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'tab':
        print("Starting full animation...")
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

# ------------------------- Show ----------------------------
print("\n--- Controls ---")
print("Press 'Enter' to advance frame.")
print("Press 'Backspace' to go back a frame.")
print("Press 'Tab' to play full animation.")
print("----------------")

animate(on_key.frame)
plt.show()
