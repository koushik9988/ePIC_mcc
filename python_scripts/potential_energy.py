import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin

file_name = 'result.h5'
path = sys.argv[1]

output_dir = './plots'
path_fig = pjoin(path, output_dir)

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')


metadata_group = f['/metadata']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
normscheme = metadata_group.attrs['norm_scheme']

mfactor = wpi / wpe
    # Get time and potential energy
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor # Time
spno = f["/metadata"].attrs["spno"]
pe = data[:, 3 * spno + 1]  # Potential energy column



plt.figure(figsize=(7, 5))
plt.plot(ts, pe, label="Potential Energy (neutral density = 5e19)", color='purple')
plt.yscale("log")
#plt.ylim([1e-4,2e0])
plt.xlabel(r"$\omega_{pi}t$")
plt.ylabel("Potential Energy")
plt.title("Potential Energy vs Time")
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()

plt.savefig(pjoin(path_fig,'pe.pdf'),dpi=1200)

plt.show()
