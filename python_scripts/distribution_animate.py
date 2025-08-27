import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# ------------------------- Argument Handling ----------------------------
if len(sys.argv) != 3:
    print("Usage: python3 distibution_animate.py <species_name>")
    sys.exit(1)

data_path = sys.argv[1] 
species_name = sys.argv[2]

# Path to the directory containing the HDF5 file
file_name = "result.h5"
path = pjoin(data_path, file_name)

# ------------------------- HDF5 Loading ----------------------------
f = h5py.File(path, 'r')
metadata = f['/metadata']

# ------------------------- Constants ----------------------------
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
normscheme = metadata.attrs['norm_scheme']
write_interval_phase = metadata.attrs['write_int_phase']
NUM_TS = metadata.attrs['NUM_TS']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

ts = f["time_var/kinetic_energy"][:, 0] * mfactor

# ------------------------- Get Species Weight ----------------------------
species_metadata = f.get(f'/metadata_species/{species_name}')

species_density = species_metadata.attrs['density']
species_spwt = species_metadata.attrs['spwt']

# ------------------------- Find Velocity Range ----------------------------
v_min, v_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    v_min = min(v_min, np.min(data_species))
    v_max = max(v_max, np.max(data_species))

# ------------------------- Initial Distribution ----------------------------
j0 = 0
data_species_0 = f[f"particle_{species_name}/{j0}"][:, 1]
nbins = 100
hist_species_0, bins = np.histogram(data_species_0, bins=nbins, range=(v_min, v_max), density=True)
hist_species_0 *= species_spwt
bin_centers_0 = 0.5 * (bins[:-1] + bins[1:])

# ------------------------- Plotting Setup ----------------------------
fig, ax = plt.subplots()

def animate(i):
    j = i * write_interval_phase
    data_species = f[f"particle_{species_name}/{j}"][:, 1]
    hist_species, bins = np.histogram(data_species, bins=nbins, range=(v_min, v_max), density=True)
    hist_species *= species_spwt
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    ax.clear()
    ax.plot(bin_centers_0, hist_species_0, linestyle='--', color='gray', label='Initial')
    ax.plot(bin_centers, hist_species, color='blue', label='Current')
    ax.set_title(f"{species_name} Distribution at $\\omega_{{pi}}t$ = {ts[i]:.2f}")
    ax.set_xlabel('Velocity (v)')
    ax.set_ylabel('f(v)')
    ax.legend()

    return ax

# ------------------------- Interaction via Keys ----------------------------
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

# ------------------------- Show ----------------------------
animate(on_key.frame)
plt.show()
