import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py

# hdf5 file name and path (assuming provided as command-line arguments)
file_name = 'result.h5'
path = sys.argv[1]
timestep = int(sys.argv[2]) if len(sys.argv) > 2 else None

# Read hdf5 file
f = h5py.File(os.path.join(path, file_name), 'r')
metadata_group = f['/metadata']

# Read attributes (assuming constants like eps0, kb, me, AMU, e are defined elsewhere)
e = constants('elementary charge')

# Constants and data loading from input.ini file
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
Te = e * metadata_group.attrs['Te']
alp = metadata_group.attrs['alpha']
n0 = metadata_group.attrs['density']

# Initialize arrays
time_steps = range(0, NUM_TS, write_interval_phase)
max_den1_values = []

# Iterate over time steps
for i in time_steps:
    den1 = f[f"fielddata/efield/{i}"][:]
    max_den1_values.append(np.max(den1))

# Plotting
fig, ax = plt.subplots()
ax.plot(time_steps, max_den1_values, color='blue', label='Max field value')
ax.set_xlabel('Time Step')
ax.set_ylabel('Max field value')
ax.legend(loc='upper right', framealpha=0.5)
plt.title('decay/growth')
plt.grid(True)
plt.show()
