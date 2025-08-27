import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os
from os.path import join as pjoin
import sys
import h5py

# ---------------- Arguments ----------------
file_name = 'result.h5'
path = sys.argv[1]  # First command-line argument

# ---------------- Plot Output Path ----------------
plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# ---------------- Open HDF5 File ----------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# ---------------- Physical Constants ----------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# ---------------- Simulation Parameters ----------------
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']


timesteps = range(0, NUM_TS + 1, write_interval)
num_timesteps = len(timesteps)


num_cells = len(f[f"fielddata/den_electron/0"]) #spatial grid size/cell no
sum_den_e = np.zeros(num_cells)
sum_den_i = np.zeros(num_cells)

# ---------------- average over all timesteps ----------------
for j in timesteps:
    den_e = f[f"fielddata/den_electron/{j}"][:]
    den_i = f[f"fielddata/den_ion/{j}"][:]
    
    sum_den_e += den_e
    sum_den_i += den_i

avg_den_e = sum_den_e / num_timesteps
avg_den_i = sum_den_i / num_timesteps


x = np.linspace(0, NC, num_cells)

plt.figure(figsize=(8, 4))
plt.plot(x, avg_den_e, label="Average Electron Density", color='black')
plt.plot(x, avg_den_i, label="Average Ion Density", color='red')
plt.xlabel("Position (x)")
plt.ylabel("Density")
plt.title("Time-Averaged Density vs Position")
plt.legend()
plt.grid(True)

plt.tight_layout()
fig_name = "avg_density_vs_position.png"
plt.savefig(pjoin(path_fig, fig_name), dpi=300)
plt.show()
