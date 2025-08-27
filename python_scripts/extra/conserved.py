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



wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']

mfactor = wpi/wpe

data = f["time_var/momentum"]
data1 = f["time_var/kinetic_energy"]

kee = data1[:,1]
kei = data1[:,2]
ken = data1[:,3]
keb = data1[:,4]
PE = data1[:,5]

ts = data[:,0]*mfactor
m_e = data[:,1]
m_i = data[:,2]
m_n = data[:,3]
m_b = data[:,4]



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

#fig, ax = plt.subplots(figsize=figsize/10.4,constrained_layout=True,dpi=ppi)
# Create a 2x3 grid of subplots
fig, ax = plt.subplots(1,2, figsize=figsize/10.4,constrained_layout=True,dpi=ppi)


# Plot on each subplot
ax[0].plot(ts, m_e, label='$P_{e}$')
ax[0].plot(ts, m_i, label='$P_{i}$')
ax[0].plot(ts, m_n, label='$P_{n}$')
ax[0].plot(ts, m_b, label='$P_{b}$')
ax[0].plot(ts, m_e + m_i + m_n + m_b, label='Total Momentum')
ax[0].set_xlabel('$\omega_{pi}t$')
ax[0].set_ylabel('$Momentum$')
#ax[0].set_ylim([min(m_e + m_i + m_n + m_b) - 10.0, max(m_e + m_i + m_n + m_b ) + 10.0])
ax[0].grid(True)
ax[0].legend(loc='upper right', framealpha=0.5)



ax[1].plot(ts, kei + kee + keb + ken, label='$KE$',color='black')
ax[1].plot(ts, PE, label="$PE$",color='red')
ax[1].plot(ts, kei + kee + keb + ken + PE, label="$KE+PE$",color='green')
#ax4.set_ylim([min(kei + kee + keb + ken + PE) - 1.0, max(kei + kee + keb + ken + PE ) + 1.0])
ax[1].set_xlabel('$\omega_{pi}t$')
ax[1].set_ylabel('$energy$')
ax[1].grid(True)
#ax4.set_ylim([min(kei + kee + keb + ken + PE) - 1.0, max(kei + kee + keb + ken + PE ) + 1.0])
ax[1].legend(loc='upper right',framealpha=0.5)

plt.tight_layout()


# Adjust layout and spacing
fig.tight_layout()


plt.savefig(pjoin(path,'energy_momentum_pe.pdf'),dpi = dpi)

plt.show()
