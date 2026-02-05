import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

file_name = 'result.h5'
path = sys.argv[1]
output_dir = './plots'
path_fig = pjoin(path, output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

#data = f["time_var/avg_collision_freq"]
data = f["time_var/electronegativity"]
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
normscheme = metadata_group.attrs['norm_scheme']

mfactor = wpi / wpe
#if normscheme == 2 or normscheme == 1:
#    mfactor = 1

time = data[:, 0]*mfactor
coll_freq = data[:, 1]

fig, ax = plt.subplots()
ax.plot(time, coll_freq, label = "neutral density = 5e19")
ax.set_xlabel(r"$\omega_{pi} t$")
#ax.set_ylabel(r"$pe$")
ax.set_ylabel(r"$\alpha$")
ax.legend()
plt.title("Time Evolution of electronagtivity")
plt.grid(True)
plt.tight_layout()

plt.savefig(pjoin(path_fig,'alpha_evolution.pdf'),dpi=1200)

plt.show()
