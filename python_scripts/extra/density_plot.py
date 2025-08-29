import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py
import matplotlib.animation as animation

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
timestep = int(sys.argv[2]) if len(sys.argv) > 2 else None

plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from input.ini file
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
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = e * metadata_group.attrs['Te']
Ti = e * metadata_group.attrs['Ti']
Tb = e * metadata_group.attrs['Tb']
Tn = e * metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Debye length calculation
ni0 = n0
ne0 = n0 / (1 + alp + beta)
ni0 = n0
nn0 = alp * ne0
nb0 = beta * ne0

LD = np.sqrt(eps0 * Te / (ne0 * e ** 2))

# Plotting
fig, ax = plt.subplots()

def animate(i):
    j = i * write_interval_phase

    den1 = f[f"fielddata/den_electron/{j}"]
    den2 = f[f"fielddata/den_ion/{j}"]
    x = np.linspace(0, NC, len(den1))
    ax.plot(x, den1[:], color='black', label="electron density")
    ax.plot(x, den2[:], color='red', label="ion density")
    ax.set_ylabel('$\phi$')
    ax.legend(loc='upper right', framealpha=0.5)

    return ax

if timestep is not None:
    # Show static plot for the specified time step
    j = timestep
    den1 = f[f"fielddata/den_electron/{j}"]
    den2 = f[f"fielddata/den_ion/{j}"]
    #den3 = f[f"fielddata/den_negion/{j}"]
    x = np.linspace(0, NC, len(den1))
    ax.plot(x, den1[:], color='black', label="electron density")
    ax.plot(x, den2[:], color='red', label="ion density")
    #ax.plot(x, den3[:], color='green', label="negative ion density")
    ax.legend(loc='upper right', framealpha=0.5)
    plt.show()
else:
    # Animate as usual
    ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=500, repeat=False)
    plt.show()
