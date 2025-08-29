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
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = e*metadata_group.attrs['Te']
Ti = e*metadata_group.attrs['Ti']
Tb = e*metadata_group.attrs['Tb']
Tn = e*metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
EV_TO_K = 11604.52 

#-----------------------Debye lenght calculation--------
ni0 = n0
ne0 = n0/(1+alp+beta)
ni0 = n0
nn0 = alp*ne0  
nb0 = beta*ne0


# Characteristic Debye length
LD = np.sqrt(eps0*Te/(ne0*e**2)) 
LDi = np.sqrt(eps0*Ti/(ni0*e**2)) 
electron_spwt = (ne0*NC*LD)/(nParticlesE)
print(electron_spwt)
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency
#---------------------------------------------------------

data = f["time_var/momentum"]



ts = data[:,0]
m_e = data[:,1]
m_i = data[:,2]
m_n = data[:,3]
m_b = data[:,4]

# Create a 2x3 grid of subplots
fig, ax = plt.subplots(2, 3, figsize=(12, 8))

# Flatten the axes array for easier plotting
ax = ax.flatten()

# Plot on each subplot
ax[0].plot(ts, m_e, label='$P_{e}$')
ax[0].set_xlabel('$\omega_{pe}t$')
ax[0].set_ylabel('$P_{e}$')
ax[0].grid(True)
ax[0].legend(loc='upper right', framealpha=0.5)

ax[1].plot(ts, m_i, label='$P_{i}$')
ax[1].set_xlabel('$\omega_{pe}t$')
ax[1].set_ylabel('$P_{i}$')
ax[1].grid(True)
ax[1].legend(loc='upper right', framealpha=0.5)

ax[2].plot(ts, m_n, label='$P_{n}$')
ax[2].set_xlabel('$\omega_{pe}t$')
ax[2].set_ylabel('$P_{n}$')
ax[2].grid(True)
ax[2].legend(loc='upper right', framealpha=0.5)

ax[3].plot(ts, m_b, label='$P_{b}$')
ax[3].set_xlabel('$\omega_{pe}t$')
ax[3].set_ylabel('$P_{b}$')
ax[3].grid(True)
ax[3].legend(loc='upper right', framealpha=0.5)

ax[4].plot(ts, m_e + m_i + m_n + m_b, label='Total Momentum')
ax[4].set_xlabel('$\omega_{pe}t$')
ax[4].set_ylabel('Momentum')
ax[4].set_ylim([min(m_e + m_i + m_n + m_b) - 10.0, max(m_e + m_i + m_n + m_b ) + 10.0])
ax[4].grid(True)
ax[4].legend(loc='upper right', framealpha=0.5)

# Adjust layout and spacing
fig.tight_layout()

if(save_fig == 1):
     plt.savefig(pjoin(path,'ke_pe.png'),dpi = 1200)

plt.show()
