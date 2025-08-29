
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# --- Command-line arguments ---
if len(sys.argv) != 4:
    print("Usage: python script.py <path> <species1> <species2>")
    sys.exit(1)

path = sys.argv[1]
species1 = sys.argv[2]
species2 = sys.argv[3]

# --- Paths and HDF5 file ---
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# --- Constants ---
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# --- Simulation parameters ---
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:", mfactor)


# --- Time axis ---
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor  # omega_pi * t

# --- Get species weights ---
metadata_s1 = f.get(f'/metadata_species/{species1}')
metadata_s2 = f.get(f'/metadata_species/{species2}')
if metadata_s1 is None or metadata_s2 is None:
    print(f"Metadata for '{species1}' or '{species2}' not found.")
    sys.exit(1)

spwt1 = metadata_s1.attrs['spwt']
spwt2 = metadata_s2.attrs['spwt']

# --- Determine velocity bounds ---
v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    v1 = f[f"particle_{species1}/{j}"][:, 1]
    v2 = f[f"particle_{species2}/{j}"][:, 1]
    combined = np.concatenate([v1, v2])
    v_min = min(v_min, np.min(combined))
    v_max = max(v_max, np.max(combined))

# --- Plot setup ---
fig, ax = plt.subplots()

# --- Initial histogram (t=0) ---
j0 = 0
v1_0 = f[f"particle_{species1}/{j0}"][:, 1]
v2_0 = f[f"particle_{species2}/{j0}"][:, 1]

nbins = 500
hist1_0, bins = np.histogram(v1_0, bins=nbins, range=(v_min, v_max), density=True)
hist2_0, _ = np.histogram(v2_0, bins=nbins, range=(v_min, v_max), density=True)

hist1_0 *= spwt1
hist2_0 *= spwt2
hist_combined_0 = (hist1_0 + hist2_0) / 2
bin_centers_0 = (bins[:-1] + bins[1:]) / 2

# --- Animation function ---
def animate(i):
    j = i * write_interval_phase
    v1 = f[f"particle_{species1}/{j}"][:, 1]
    v2 = f[f"particle_{species2}/{j}"][:, 1]

    nbins = 100
    hist1, bins = np.histogram(v1, bins=nbins, range=(v_min, v_max), density=True)
    hist2, _ = np.histogram(v2, bins=nbins, range=(v_min, v_max), density=True)

    hist1 *= spwt1
    hist2 *= spwt2
    hist_combined = (hist1 + hist2) / 2
    bin_centers = (bins[:-1] + bins[1:]) / 2

    ax.clear()
    ax.plot(bin_centers_0, hist_combined_0, linestyle='--', color='gray', label=f"Initial Combined")
    ax.plot(bin_centers, hist_combined, color='blue', label=f"{species1} + {species2}")
    ax.set_xlabel('Velocity (v)')
    ax.set_ylabel('f(v)')
    ax.legend(loc='upper right', framealpha=0.5)
    ax.set_title(f"$\omega_{{pi}}t$ = {ts[i]:.2f}")

    return ax

# --- Keyboard controls ---
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
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

# --- Start ---
animate(on_key.frame)
plt.show()
