import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin

file_name = 'result.h5'
path1 = sys.argv[1]
path2 = sys.argv[2]
path3 = sys.argv[3]
path4 = sys.argv[4]

output_dir = './plots'
path_fig = pjoin(path1, output_dir)
os.makedirs(path_fig, exist_ok=True)

# Open HDF5 from first path to get time
f = h5py.File(pjoin(path1, file_name), 'r')
metadata_group = f['/metadata']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
mfactor = wpi / wpe
data = f["time_var/kinetic_energy"]
data1 = f["time_var/electronegativity"]
ts = data[:, 0] * mfactor  # Time
spno = f["/metadata"].attrs["spno"]
f.close()

def load_pe(path):
    f = h5py.File(pjoin(path, file_name), 'r')
    data = f["time_var/kinetic_energy"]
    data1 = f["time_var/electronegativity"]
    spno = f["/metadata"].attrs["spno"]
    pe = data[:, 3 * spno + 1]  # Potential energy column
    #coll_freq = data1[:, 1]
    f.close()
    return pe
    #return coll_freq

# Load potential energy from all 4 paths
pe1 = load_pe(path1)
pe2 = load_pe(path2)
pe3 = load_pe(path3)
pe4 = load_pe(path4)



# Plot together
fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(ts, pe1, label="Neutral density = 1e18", color='purple')
ax.plot(ts, pe2, label="Neutral density = 5e18", color='blue')
ax.plot(ts, pe3, label="Neutral density = 1e19", color='green')
ax.plot(ts, pe4, label="Neutral density = 5e19", color='red')

#ax.set_yscale("log")
#ax.set_ylim([1e-4, 2e0])
ax.set_xlabel(r"$\omega_{pi} t$")
ax.set_ylabel("Potential Energy")
#ax.set_ylabel(r"$\alpha$")
#ax.set_title("electronegativity over time")
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend()
plt.tight_layout()

plt.savefig(pjoin(path_fig, 'pe_comp.pdf'), dpi=1200)
plt.show()
