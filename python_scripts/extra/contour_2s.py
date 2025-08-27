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
        phi = f['fielddata/pot/' + str(time_step)]
        x = np.linspace(0,NC,len(phi))
        pot[i, :] = phi
        EF  = f['fielddata/efield/' + str(time_step)]
        efield[i,:] = EF
        nde = f['fielddata/den_electron/' + str(time_step)]
        eden[i,:] = nde
        ndi = f['fielddata/den_ion/' + str(time_step)]
        iden[i,:] = ndi



phi_cbrt = np.sign(pot) * np.abs(pot)**(1/3)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(15, 8))
ax1.contourf(x, ts, phi_cbrt)
ax1.set_xlabel('x')
ax1.set_ylabel('$\omega_{pi}t$')
ax1.set_title('Electric Potential')

ax2.contourf(x, ts, efield)
ax2.set_xlabel('x')
ax2.set_ylabel('$\omega_{pi}t$')
ax2.set_title('Electric Field')

ax3.contourf(x, ts, eden)
ax3.set_xlabel('x')
ax3.set_ylabel('$\omega_{pi}t$')
ax3.set_title('Electron Density')

ax4.contourf(x, ts, iden)
ax4.set_xlabel('x')
ax4.set_ylabel('$\omega_{pi}t$')
ax4.set_title('Ion Density')


plt.show()
