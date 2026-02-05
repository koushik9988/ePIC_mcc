import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

# HDF5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
timestep = int(sys.argv[2]) if len(sys.argv) > 2 else None

plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Prepare data for contour plot
x = np.linspace(0, NC, NC)  # Spatial grid
t = np.arange(0, DATA_TS_PHASE) * write_interval_phase * DT_coeff  # Time steps
temperature = np.zeros((DATA_TS_PHASE, NC))  # 2D array for temperature

# Compute temperature-like quantity (v^2 as proxy for temperature)
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    vel_electron = f[f"fielddata/vel_electron/{j}"][:]  # Electron velocity
    temperature[i, :] = vel_electron**2  # Proportional to kinetic energy

# Create contour plot
fig, ax = plt.subplots()
X, T = np.meshgrid(x, t)  # Create 2D grid for contour
contour = ax.contourf(X, T, temperature, cmap='hot', levels=50)
fig.colorbar(contour, label='Temperature (proportional to $v_e^2$)')
ax.set_xlabel('x')
ax.set_ylabel('Time (s)')
ax.set_title('Temperature Contour Plot')

plt.savefig(pjoin(path_fig, 'temperature_contour.png'))
plt.show()