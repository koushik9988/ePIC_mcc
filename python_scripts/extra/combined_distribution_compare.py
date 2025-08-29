import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Set paths for two data files
path1 = sys.argv[1]#"../data_coll_1e15" 
path2 = sys.argv[2]#"../data_nocoll"  

# HDF5 file name and paths for plots
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path1, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 files
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes (assuming they are the same for both files; adjust if different)
NC = metadata_group1.attrs['NC']
NUM_TS = metadata_group1.attrs['NUM_TS']
write_interval = metadata_group1.attrs['write_int']
write_interval_phase = metadata_group1.attrs['write_int_phase']
DT_coeff = metadata_group1.attrs['DT_coeff']
LD = metadata_group1.attrs['LDe']
LDi = metadata_group1.attrs['LDi']
wpe = metadata_group1.attrs['wpe']
wpi = metadata_group1.attrs['wpi']
save_fig = metadata_group1.attrs['save_fig']
normscheme = metadata_group1.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:", mfactor)

data1 = f1["time_var/kinetic_energy"]
ts = data1[:, 0]
ts *= mfactor  # converted to w_pi*t

# Get specific weights for electron and beam species from both files
metadata_electron1 = f1.get('/metadata_species/electron')
metadata_beam1 = f1.get('/metadata_species/beam')
metadata_electron2 = f2.get('/metadata_species/electron')
metadata_beam2 = f2.get('/metadata_species/beam')

if not (metadata_electron1 and metadata_beam1 and metadata_electron2 and metadata_beam2):
    print("Metadata for 'electron' or 'beam' not found in one or both files.")
    sys.exit(1)

electron_spwt1 = metadata_electron1.attrs['spwt']
beam_spwt1 = metadata_beam1.attrs['spwt']
electron_spwt2 = metadata_electron2.attrs['spwt']
beam_spwt2 = metadata_beam2.attrs['spwt']

# Find velocity range across both species and both files
v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_electron1 = f1[f"particle_electron/{j}"][:, 1]
    data_beam1 = f1[f"particle_beam/{j}"][:, 1]
    data_electron2 = f2[f"particle_electron/{j}"][:, 1]
    data_beam2 = f2[f"particle_beam/{j}"][:, 1]
    combined_data = np.concatenate([data_electron1, data_beam1, data_electron2, data_beam2])
    v_min = min(v_min, np.min(combined_data))
    v_max = max(v_max, np.max(combined_data))

# Plotting setup
fig, ax = plt.subplots()

# Compute initial histogram (t=0) for combined data from both files
j0 = 0  # Initial time step
data_electron1_0 = f1[f"particle_electron/{j0}"][:, 1]
data_beam1_0 = f1[f"particle_beam/{j0}"][:, 1]
data_electron2_0 = f2[f"particle_electron/{j0}"][:, 1]
data_beam2_0 = f2[f"particle_beam/{j0}"][:, 1]

nbins = 500
hist_electron1_0, bins = np.histogram(data_electron1_0, bins=nbins, range=(v_min, v_max), density=True)
hist_beam1_0, _ = np.histogram(data_beam1_0, bins=nbins, range=(v_min, v_max), density=True)
hist_electron2_0, _ = np.histogram(data_electron2_0, bins=nbins, range=(v_min, v_max), density=True)
hist_beam2_0, _ = np.histogram(data_beam2_0, bins=nbins, range=(v_min, v_max), density=True)

# Scale each histogram by its specific weight
hist_electron1_0 *= electron_spwt1
hist_beam1_0 *= beam_spwt1
hist_electron2_0 *= electron_spwt2
hist_beam2_0 *= beam_spwt2

# Combine the scaled histograms for each file
hist_combined1_0 = (hist_electron1_0 + hist_beam1_0) / 2
hist_combined2_0 = (hist_electron2_0 + hist_beam2_0) / 2
bin_centers_0 = (bins[:-1] + bins[1:]) / 2

def animate(i):
    j = i * write_interval_phase
    data_electron1 = f1[f"particle_electron/{j}"][:, 1]
    data_beam1 = f1[f"particle_beam/{j}"][:, 1]
    data_electron2 = f2[f"particle_electron/{j}"][:, 1]
    data_beam2 = f2[f"particle_beam/{j}"][:, 1]

    # Histograms for species velocities with separate scaling
    nbins = 100
    hist_electron1, bins = np.histogram(data_electron1, bins=nbins, range=(v_min, v_max), density=True)
    hist_beam1, _ = np.histogram(data_beam1, bins=nbins, range=(v_min, v_max), density=True)
    hist_electron2, _ = np.histogram(data_electron2, bins=nbins, range=(v_min, v_max), density=True)
    hist_beam2, _ = np.histogram(data_beam2, bins=nbins, range=(v_min, v_max), density=True)

    # Scale each histogram by its specific weight
    hist_electron1 *= electron_spwt1
    hist_beam1 *= beam_spwt1
    hist_electron2 *= electron_spwt2
    hist_beam2 *= beam_spwt2

    # Combine the scaled histograms for each file
    hist_combined1 = (hist_electron1 + hist_beam1) / 2
    hist_combined2 = (hist_electron2 + hist_beam2) / 2
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Clear and plot
    ax.clear()

    # Plot initial combined distributions with dashed lines
    ax.plot(bin_centers_0, hist_combined1_0, linestyle='--', color='gray', 
            label=f"Initial Combined ({path1})")
    ax.plot(bin_centers_0, hist_combined2_0, linestyle='--', color='black', 
            label=f"Initial Combined ({path2})")

    # Plot combined distributions at current time step
    ax.plot(bin_centers, hist_combined1, color='blue', 
            label=f"Combined (Electron + Beam) ({path1})")
    ax.plot(bin_centers, hist_combined2, color='red', 
            label=f"Combined (Electron + Beam) ({path2})")
    
    ax.set_xlabel('Velocity (v)')
    ax.set_ylabel('f(v)')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_title(f'$\omega_{{pi}}t$ = {ts[i]:.2f} (Combined Distribution Comparison)')

    return ax

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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, 
                                             blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None  # To keep reference to animation

fig.canvas.mpl_connect('key_press_event', on_key)

# Start with the first frame
animate(on_key.frame)
plt.show()