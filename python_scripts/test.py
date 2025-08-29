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
ts = data[:, 0] * mfactor  # Time
f.close()

def load_pe(path):
    f = h5py.File(pjoin(path, file_name), 'r')
    gas_density = f['/metadata'].attrs['GAS_DENSITY']
    data = f["time_var/kinetic_energy"]
    spno = f["/metadata"].attrs["spno"]
    pe = data[:, 3 * spno + 1]  # Potential energy column
    f.close()
    return pe, gas_density

# Load potential energy + density from all 4 paths
pe1, gd1 = load_pe(path1)
pe2, gd2 = load_pe(path2)
pe3, gd3 = load_pe(path3)
pe4, gd4 = load_pe(path4)

figsize = np.array([80,80/1.618])
ppi = np.sqrt(1920**2+1200**2)/24
fig, ax = plt.subplots(figsize=figsize/10.4, constrained_layout=True, dpi=ppi)

ax.plot(ts, pe1, label=fr"Neutral density = {gd1:.2e}", color='purple')
ax.plot(ts, pe2, label=fr"Neutral density = {gd2:.2e}", color='blue')
ax.plot(ts, pe3, label=fr"Neutral density = {gd3:.2e}", color='green')
ax.plot(ts, pe4, label=fr"Neutral density = {gd4:.2e}", color='red')

ax.set_xlabel(r"$\omega_{pi} t$")
ax.set_ylabel("Potential Energy")
ax.semilogy()
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend()
plt.tight_layout()

plt.savefig(pjoin(path_fig, 'pe_comp.pdf'), dpi=1200)
plt.show()
