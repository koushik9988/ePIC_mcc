import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py
from scipy.signal import find_peaks

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
NC = metadata_group.attrs['NC']
save_fig = metadata_group.attrs['save_fig']
fielddata_group = f['fielddata/pot']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
n0 = metadata_group.attrs['density']
write_interval = metadata_group.attrs['write_int']
NUM_TS = metadata_group.attrs['NUM_TS']
ne0 = n0
we = np.sqrt(ne0 * e**2 / (eps0 * me))

# Fixed spatial point
fixed_point_index = 1

data = f["time_var/kinetic_energy"]
ts = data[:, 0]
times = ts * (1 / we)

num_timesteps = int(NUM_TS / write_interval)
print(num_timesteps)
potentials = np.zeros(num_timesteps+1)

for i in range(num_timesteps+1):
    j = i * write_interval
    #print(i)
    pot = f[f"fielddata/pot/{j}"]
    #pot = f[f"fielddata/den_electron/{j}"]
    #pot = f[f"fielddata/den_ion/{j}"]
    potentials[i] = pot[fixed_point_index]

# Find peaks
peaks, _ = find_peaks(potentials)

peak_times = times[peaks]
peak_diffs = np.diff(peak_times)

fig, ax = plt.subplots()
ax.plot(times, potentials, color='black', label="Potential at fixed point")
ax.set_xlabel("Time")
ax.set_ylabel("Potential")
ax.legend()

# Print peak differences
print("Differences between consecutive peak positions:", peak_diffs)

plt.savefig(pjoin(path, 'ke_pe.png'), dpi=1200)
plt.show()

# Close the file
f.close()
