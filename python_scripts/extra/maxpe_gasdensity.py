import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import h5py
import os
from os.path import join as pjoin

# HDF5 file name and paths
file_name = 'result.h5'
path1 = "../data_esw_nocoll"
path2 = "../data_esw_coll_1e11"
path3 = "../data_esw_coll_1e15"
path4 = "../data_esw_coll_1e17"
path5 = "../data_esw_coll_1e19"
path6 = "../data_esw_coll_1e20"

path_fig = './plots'

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 files
paths = [path1, path2, path3, path4, path5, path6]
colors = ['red', 'green', 'blue', 'purple', 'violet', 'black']

# Manually defining gas density values for the x-axis
gas_densities = [0, 1e11, 1e15, 1e17, 1e19, 1e20]  # Corresponding to each dataset

max_PE_values = []

for path in paths:
    with h5py.File(pjoin(path, file_name), 'r') as f:
        data = f["time_var/kinetic_energy"]
        PE = data[:, 4]  # Electrostatic potential energy
        max_PE = np.max(PE)
        max_PE_values.append(max_PE)

# Plotting Max PE vs Gas Density
fig, ax = plt.subplots(figsize=(8, 6), dpi=120)
ax.plot(gas_densities, max_PE_values)
ax.set_xlabel('Gas Density')
ax.set_ylabel('Max Electrostatic Potential Energy')
#ax.grid(True, which='both', linestyle='--')
plt.savefig(pjoin(path_fig, 'max_pe_vs_gas_density.png'), dpi=300)
plt.show()