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
import scipy.integrate as intg



file_name = 'result.h5'
path = sys.argv[1]
output_dir = './plots'
path_fig = pjoin(path, output_dir)
#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']


# Constants and data loading from hdf5 file 
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read individual attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
EV_TO_K = 11604.52 

mFactor = wpi/wpe
DATA_TS = int(NUM_TS/write_interval) + 1
# ------------------------------------------------------
L = LD
data = f["time_var/kinetic_energy"]
ts = data[:,0]
ts *= mFactor # This is time w_{pi}t
#------------- potential-energy calculation-------------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
PE = np.zeros(len(time_steps))

pot = np.empty(shape=(DATA_TS, NC))
efield = np.empty(shape=(DATA_TS, NC))
eden = np.empty(shape=(DATA_TS, NC))
iden = np.empty(shape=(DATA_TS, NC))

# ------------------------------------------------------
for i, time_step in enumerate(time_steps):
        phi = f['fielddata/coll_rate_ion/' + str(time_step)]
        x = np.linspace(0,NC,len(phi))
        pot[i, :] = phi
        EF  = f['fielddata/efield/' + str(time_step)]
        efield[i,:] = EF
       


fig, ax = plt.subplots(figsize=(12, 6))
X, Y = np.meshgrid(np.linspace(0, NC, pot.shape[1]), ts)
cf = ax.contourf(X, Y, pot, levels=50, cmap='plasma')

ax.set_xlabel('x')
ax.set_ylabel(r'$\omega_{pi} t$')
ax.set_title('coll_rate (filled contours)')
fig.colorbar(cf, ax=ax, label='coll_rate')
plt.show()


"""
fig, ax = plt.subplots(figsize=(12, 6))

im = ax.imshow(pot,aspect='auto', origin='lower',extent=[0, NC, ts[0], ts[-1]])

ax.set_xlabel('x')
ax.set_ylabel(r'$\omega_{pi} t$')
ax.set_title('coll_rate')
fig.colorbar(im, ax=ax, label='data')

plt.savefig(pjoin(path_fig,'coll_rate_negion.pdf'),dpi=1200)
plt.show()
"""