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


#hdf5 file name and path 
file_name = 'result.h5'
path = sys.argv[1]

path1 = './plots'

path_fig = pjoin(path,path1)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

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

data = f["time_var/kinetic_energy"]

mFactor = wpi/wpe

ts = data[:,0]
kee = data[:,1]
kei = data[:,2]
#ken = data[:,3]
PE = data[:,3]


if normscheme == 2 or normscheme == 4:
    LD = LDi
    ts *= mFactor
if normscheme == 1 or normscheme == 3:
    LD = LD

#---------------plotting -------------------------
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

fig, ((ax1, ax2),(ax3, ax4))= plt.subplots(2, 2, figsize= figsize/10.4,constrained_layout=True,dpi=ppi)

ax1.plot(ts, kee, label='$KE_{e}$')
ax1.set_xlabel('$\omega_{pi}t$')
ax1.set_ylabel('$KE_{e}$')
ax1.grid(True)
ax1.legend(loc='upper right',framealpha=0.5)

ax2.plot(ts, kei, label="$KE_{i}$")
ax2.set_xlabel('$\omega_{pi}t$')
ax2.set_ylabel('$KE_{i}$')
ax2.grid(True)
ax2.legend(loc='upper right',framealpha=0.5)

#ax5.plot(ts, ken, label="$KE_{n}$")
#ax5.set_xlabel('$\omega_{pi}t$')
#ax5.set_ylabel('$KE_{n}$')
#ax5.grid(True)
#ax5.legend(loc='upper right',framealpha=0.5)



ax3.plot(ts, PE, label="potential energy")
ax3.set_xlabel('$\omega_{pi}t$')
ax3.set_ylabel('$PE$')
ax3.grid(True)
ax3.legend(loc='upper right',framealpha=0.5)

ax4.plot(ts, kei + kee + PE, label="Total Energy")
ax4.set_xlabel('$\omega_{pi}t$')
ax4.set_ylabel('Total Energy')
ax4.grid(True)
ax4.set_ylim([min(kei + kee  + PE) - 1, max(kei + kee  + PE ) + 1])
ax4.legend(loc='upper right',framealpha=0.5)

plt.tight_layout()

#if(save_fig == 1):
plt.savefig(pjoin(path,'ke_pe.png'),dpi = dpi)

plt.show()