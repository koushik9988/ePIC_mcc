import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Set path
path = sys.argv[1]  # "../data_esw_coll_1e19"

# HDF5 file name and paths for plots
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval_phase = metadata_group.attrs['write_int_phase']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
normscheme = metadata_group.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1

# Time series data
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor  # converted to w_pi*t

# Get species metadata
metadata_electron = f.get('/metadata_species/electron')
metadata_beam = f.get('/metadata_species/beam')

electron_spwt = metadata_electron.attrs['spwt']
beam_spwt = metadata_beam.attrs['spwt']

# Find velocity range
v_min, v_max = float('inf'), float('-inf')
x_min, x_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_electron = f[f"particle_electron/{j}"]
    data_beam = f[f"particle_beam/{j}"]
    combined_data = np.concatenate([data_electron[:, 1], data_beam[:, 1]])
    v_min, v_max = min(v_min, np.min(combined_data)), max(v_max, np.max(combined_data))
    x_min, x_max = min(x_min, np.min(data_electron[:, 0])), max(x_max, np.max(data_electron[:, 0]))

# Plot setup
fig, (ax_phase, ax_dist) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

def animate(i):
    j = i * write_interval_phase
    data_electron = f[f"particle_electron/{j}"]
    data_beam = f[f"particle_beam/{j}"]

    # Phase space plot (velocity vs. position)
    ax_phase.clear()
    ax_phase.scatter(data_electron[:, 1], data_electron[:, 0], s=1, color='blue', alpha=0.5, label='Electrons')
    ax_phase.scatter(data_beam[:, 1], data_beam[:, 0], s=1, color='red', alpha=0.5, label='Beam')
    ax_phase.set_xlabel('Velocity (v)')
    ax_phase.set_ylabel('Position (x)')
    ax_phase.set_xlim(v_min, v_max)
    ax_phase.set_ylim(x_min, x_max)
    ax_phase.legend()
    ax_phase.set_title(f'Phase Space ($\omega_{{pi}}t$ = {ts[i]:.2f})')

    # Velocity distribution function
    nbins = 500
    hist_electron, bins = np.histogram(data_electron[:, 1], bins=nbins, range=(v_min, v_max), density=True)
    hist_beam, _ = np.histogram(data_beam[:, 1], bins=nbins, range=(v_min, v_max), density=True)
    hist_electron *= electron_spwt
    hist_beam *= beam_spwt
    hist_combined = (hist_electron + hist_beam) / 2
    bin_centers = (bins[:-1] + bins[1:]) / 2

    ax_dist.clear()
    ax_dist.plot(bin_centers, hist_combined, color='blue', label="Combined Distribution")
    ax_dist.set_xlabel('Velocity (v)')
    ax_dist.set_ylabel('f(v)')
    ax_dist.legend()
    ax_dist.set_title('Velocity Distribution')
    
    return ax_phase, ax_dist

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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None  # To keep reference to animation

fig.canvas.mpl_connect('key_press_event', on_key)

# Start with the first frame
animate(on_key.frame)
plt.show()

"""
fig.canvas.mpl_connect('key_press_event', lambda event: animate(min(max(0, int(event.key) if event.key.isdigit() else 0), DATA_TS_PHASE - 1)))
ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, interval=100, repeat=False)
plt.show()
"""