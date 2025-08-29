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
n = sys.argv[2]
n = int(n)
#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_negion = f['/metadata_species/negion']
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
nParticlesE = 500000#metadata_group.attrs['nE']
nParticlesI = 500000#metadata_group.attrs['nI']
nParticlesN = 500000#metadata_group.attrs['nN']
nParticlesB = 500000#metadata_group.attrs['nB']
Te = 1*e
Ti = 0.026*e
Tb = 0.026*e
Tn = 0.026*e
alp = metadata_negion.attrs['density']
beta = 0.4
mi = 1#metadata_group.attrs['mI']
mn = 1#metadata_group.attrs['mN']
mb = 1#metadata_group.attrs['mB']
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
        #x = f['fielddata/position/' + str(time_step)]
        #dx = x[1]-x[0]
        #xl = np.max(x)

xl = NC
electron_spwt = (ne0*xl*L)/(nParticlesE)
ion_spwt = (ni0*xl*L)/(nParticlesI)
neg_spwt = (nn0*xl*L)/(nParticlesN)
beam_spwt = (nb0*xl*L)/(nParticlesB)

# ========================================================================
# File index value is k(say), change this index to get the plot at required time.
# Calculate the wpet for a k-value by using DT_coeff beforehand.
k = [0,n]
nbins =500#0

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
    data_phase_e = f["particle_electron/%d"%j]
    #data_phase_e = f["particle_beam/%d"%j]  
    xe = data_phase_e[:,0]
    ve = data_phase_e[:,1]

    # Load the ion data
    data_phase_i = f["particle_ion/%d"%j]    
    xi = data_phase_i[:,0]
    vi = data_phase_i[:,1]

    # Load the negative ion data
    data_phase_n = f["particle_negion/%d"%j]    
    xn = data_phase_n[:,0]
    vn = data_phase_n[:,1]

    # Load the beam data
    data_phase_b = f["particle_beam/%d"%j]    
    xb = data_phase_b[:,0]
    vb = data_phase_b[:,1]
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

    v = np.concatenate((vn, vb))
    # +++++++++++++++++++++++++++++++++++++++++++++++++++   
    #v = np.concatenate((vi, vn, vb)) 
    # v = np.concatenate((vi, vb))

    mini = np.min(v)-0.5
    maxi = np.max(v)+0.5
    #nbins = 51
    delta = (maxi - mini)/(nbins)
    bin = np.zeros(nbins)

    # If there are three species: here 100000 is the no of particles
    for j in range(len(v)):
        c = int(np.floor((v[j] - mini)/delta))
        if j<500000:    
            bin[c] += 1*neg_spwt
        else:
            bin[c] += 1*beam_spwt
    
    # If there are two species
    # for j in range(len(v)):
    #     c = int(np.floor((v[j] - mini)/delta))
    #     if j<100000:    
    #         bin[c] += 1*ion_spwt        
    #     else:
    #         bin[c] += 1*beam_spwt

    y = np.linspace(np.min(v), np.max(v), len(bin))    
    count[i,:] = bin
    veli[i,:] = y
    
# -----------------------------------------------------
# ------------------------------------------------------
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


fig,(ax,ax1) = plt.subplots(2,1,figsize=figsize/10.4,constrained_layout=True,dpi=ppi)
#ax.plot(veli[0,:], count[0,:]/np.max(count[0,:]),'r-', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[0]))
#ax.plot(veli[1,:], count[1,:]/np.max(count[1,:]),'k--', linewidth=2,label='$\u03C9_{pi}t=%d$'%(wpit[1]))

ax.plot(veli[0,:], count[0,:]/np.max(count[0,:]),'r-', linewidth=2,label="time_step = %d"%(n))
ax.plot(veli[1,:], count[1,:]/np.max(count[1,:]),'k--', linewidth=2,label="time_step =%d"%(n))

ax.set_xlabel('v')
ax.set_ylabel('f(v)')
ax.set_title('Combined Ion Velocity Distribution')
ax.legend(loc=0)

ax1.plot(vele[0,:], count_e[0,:]/np.max(count_e[0,:]),'r-', linewidth=2,label="time_step =%d"%(n))
ax1.plot(vele[1,:], count_e[1,:]/np.max(count_e[1,:]),'k--', linewidth=2,label="time_step =%d"%(n))
ax1.set_xlabel('v')
ax1.set_ylabel('f(v)')
ax1.set_title('Electron Velocity Distribution')
ax1.legend(loc=0)

# Save Figure
plt.savefig(pjoin(path,'veldist_spline_%d.png'%(10)),format='pdf',dpi=dpi)
plt.savefig(pjoin(path,'veldist_%d.png'%(10)),format='pdf',dpi=dpi)

plt.show()
