import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
import seaborn as sns
from scipy.constants import value as constants 
from scipy.interpolate import make_interp_spline
import h5py 

sns.set(style='whitegrid')
import os.path
from os.path import join as pjoin
import sys

#-------------------------------------------------------------------------------
file_name = 'result.h5'
path = sys.argv[1]
#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
# -----------------------------------------------------------------
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
nParticlesN = metadata_group.attrs['nN']
nParticlesB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']*e
Ti = metadata_group.attrs['Ti']*e
Tb = metadata_group.attrs['Tb']*e
Tn = metadata_group.attrs['Tn']*e
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
#save_fig = metadata_group.attrs['save_fig']
EV_TO_K = 11600 
DATA_TS = int(NUM_TS/write_interval) + 1
# ------------------------------------------------------
ni0 = n0
ne0 = n0/(1+alp+beta)
nn0 = alp*ne0  
nb0 = beta*ne0

LDE = np.sqrt(eps0*Te/(ne0*e**2)) 
LDI = np.sqrt(eps0*Ti/(ni0*e**2))

wpe = np.sqrt(ne0*e**2/(eps0*me)) # electron Plasma Frequency
wpi = np.sqrt(ni0*e**2/(eps0*mi)) # ion Plasma Frequency

mFactor = wpi/wpe # w_pi*t = x * w_pe * t => x = w_pi/w_pe 
# print(mFactor)

L = LDI
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        x = np.linspace(0,NC,len(EF_data))*L
        dx = x[1]-x[0]
        xl = np.max(x)

electron_spwt = (ne0*xl*L)/(nParticlesE)
ion_spwt = (ni0*xl*L)/(nParticlesI)

# ========================================================================
# File index value is k(say), change this index to get the plot at required time.
# Calculate the wpet for a k-value by using DT_coeff beforehand.
k = [0,50000]
nbins = 500 #50

count_e = np.zeros([len(k), nbins]) 
count = np.zeros([len(k), nbins])

vele = np.zeros([len(k), nbins])
veli = np.zeros([len(k), nbins])

wpet = np.zeros(len(k))
wpit = np.zeros(len(k))

for i in range(len(k)):   
    j = k[i]
    wpet[i] = k[i]*DT_coeff # Value of wpet
    wpit[i] = mFactor*wpet[i] # Value of wpit converted from wpet

    # Load the electron data
    data_phase_e = f["particle_electron/%d"%j ]   
    xe = data_phase_e[:,0]
    ve = data_phase_e[:,1]

    # Load the ion data
    data_phase_i = f["particle_ion/%d"%j ]    
    xi = data_phase_i[:,0]
    vi = data_phase_i[:,1]

    # +++++++++++++++++++++++++++++++++++++++++++++++++++
    # Calculations for electrons
    min_e = np.min(ve)-0.5
    max_e = np.max(ve)+0.5
    delta_e = (max_e - min_e)/(nbins)
    bin_e = np.zeros(nbins)
    for j in range(len(ve)):
        ce = int(np.floor((ve[j] - min_e)/delta_e))
        # print(ce)
        bin_e[ce] += 1*electron_spwt
    
    ye = np.linspace(np.min(ve), np.max(ve), len(bin_e))    
    count_e[i,:] = bin_e
    vele[i,:] = ye
    # +++++++++++++++++++++++++++++++++++++++++++++++++++

    # Calculations for ion
    min_i = np.min(vi)-0.5
    max_i = np.max(vi)+0.5
    delta_i = (max_i - min_i)/(nbins)
    bin_i = np.zeros(nbins)
    for j in range(len(vi)):
        ci = int(np.floor((vi[j] - min_i)/delta_i))
        # print(ce)
        bin_i[ci] += 1*ion_spwt
    
    yi = np.linspace(np.min(vi), np.max(vi), len(bin_i))    
    count_e[i,:] = bin_i
    veli[i,:] = yi
    # +++++++++++++++++++++++++++++++++++++++++++++++++++      
    

# ------------------------------------------------------
fig, ax = plt.subplots()
ax.plot(veli[0,:], count[0,:]/np.max(count[0,:]),'r-', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[0]))
ax.plot(veli[1,:], count[1,:]/np.max(count[1,:]),'k--', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[1]))
ax.set_xlabel('v')
ax.set_ylabel('f(v)')
ax.set_title('Ion Velocity Distribution')
ax.legend(loc=0)

fig, ax1 = plt.subplots()
ax1.plot(vele[0,:], count_e[0,:]/np.max(count_e[0,:]),'r-', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[0]))
ax1.plot(vele[1,:], count_e[1,:]/np.max(count_e[1,:]),'k--', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[1]))
ax1.set_xlabel('v')
ax1.set_ylabel('f(v)')
ax1.set_title('Electron Velocity Distribution')
ax1.legend(loc=0)

# Save Figure
# plt.savefig(pjoin(path_fig,'veldist_spline_%d.png'%(vd)),format='png',dpi=600)
# plt.savefig(pjoin(path_fig,'veldist_%d.png'%(vd)),format='png',dpi=600)

plt.show()
