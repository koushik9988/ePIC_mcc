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


eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes from hdf5 file 
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
#write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Ti = metadata_group.attrs['Ti']
Tb = metadata_group.attrs['Tb']
Tn = metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
#save_fig = metadata_group.attrs['save_fig']
EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1
# ------------------------------------------------------
ni0 = n0
ne0 = n0/(1+alp+beta)
nn0 = alp*ne0  
nb0 = beta*ne0
# ------------------------------------------------------
LDE = np.sqrt(eps0*Te/(ne0*e**2)) 
LDI = np.sqrt(eps0*Ti/(ni0*e**2))
wpe = np.sqrt(ne0*e**2/(eps0*me)) # electron Plasma Frequency
wpi = np.sqrt(ni0*e**2/(eps0*mi)) # ion Plasma Frequency
# ------------------------------------------------------
L = LDI
# ------------------------------------------------------
mFactor = wpi/wpe
# ------------------------------------------------------
data = f["time_var/kinetic_energy"]
ts = data[:,0]
ts *= mFactor # This is time w_{pi}t
#------------- potential-energy calculation-------------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
PE = np.zeros(len(time_steps))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        x = np.linspace(0,NC,len(EF_data))*L
        # --------------------------------------------
        dx = x[1]-x[0]
        xl = NC*dx
        # --------------------------------------------
        #E_sq = (EF_data[:]** 2) * ((me*we**2*LD/e)**2) 
        E_sq = (EF_data[:]** 2)*((n0*e*L/eps0)**2)
        integral = intg.trapz(E_sq, x)
        PE[i] = 0.5 * eps0 * integral

# electron_spwt = (ne0*NC*LDE)/(nParticlesE)
# ion_spwt = (ni0*xl*L)/(nParticlesI)

# TH = (ion_spwt*nParticlesI)*Ti 
# # Normalize the electric potential energy by the total ion thermal energy 
# PE/= TH

pot = np.empty(shape=(DATA_TS, NC+1))
efield = np.empty(shape=(DATA_TS, NC+1))
eden = np.empty(shape=(DATA_TS, NC+1))
iden = np.empty(shape=(DATA_TS, NC+1))
nden = np.empty(shape=(DATA_TS, NC+1))
bden = np.empty(shape=(DATA_TS, NC+1))
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
        ndn = f['fielddata/den_negion/' + str(time_step)]
        nden[i,:] = ndn
        ndb = f['fielddata/den_beam/' + str(time_step)]
        bden[i,:] = ndb

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(15, 8))
ax1.contourf(x, ts, pot)
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

ax5.contourf(x, ts, nden)
ax5.set_xlabel('x')
ax5.set_ylabel('$\omega_{pi}t$')
ax5.set_title('Negative Ion Density')

ax6.contourf(x, ts, bden)
ax6.set_xlabel('x')
ax6.set_ylabel('$\omega_{pi}t$')
ax6.set_title('Beam Density')

plt.savefig(pjoin(path,'contour.png'),dpi= 1200)
plt.show()
