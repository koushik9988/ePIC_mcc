"""
Plots the dispersion graph of from the electric field data
Use %reset -f to clear the python workspace
Data File Invoked: processed_results_all.npz
Run as: 
"""
# TO RUN: python3 dispersion_npz.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

file_name = 'result.h5'
#path = '../data_run_1/'
path = sys.argv[1]
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']
metadata_negion = f['/metadata_species/negion']
metadata_beam = f['/metadata_species/beam']
metadata_group = f['/metadata']

#path = sys.argv[1]
# ------------------ Comments -------------------------------------------------
# input parameters specific file
# path to the data folder is ../data/data002_vd_20/files for vd=20
# path to the data folder is ../data/data001_vd_80/files for vd=80
#------------------------------------------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read attributes from hdf5 file 

NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']

Te = metadata_electron.attrs['temperature']
Ti = metadata_ion.attrs['temperature']
Tn = metadata_negion.attrs['temperature']
Tb = metadata_beam.attrs['temperature']
mi = metadata_ion.attrs['mass']
mn = metadata_negion.attrs['mass']
mb = metadata_beam.attrs['mass']
vb = metadata_beam.attrs['streaming_velocity']
nParticlesE = metadata_electron.attrs['num_particles']
nParticlesI = metadata_ion.attrs['num_particles']
nParticlesN = metadata_negion.attrs['num_particles']
nParticlesB = metadata_beam.attrs['num_particles']
alp = metadata_negion.attrs['density']
beta = metadata_beam.attrs['density']

EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1

we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']

# print(normscheme)
#------------------------------------------------------
vthi = np.sqrt(Ti/mi)
vthe = np.sqrt(Te/me)
vthn = np.sqrt(Tn/mn)
#-------------------------------------------------------
ni0 = n0
ne0 = n0/(1+alp+beta)
nn0 = alp*ne0
nb0 = beta*ne0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

LDE = np.sqrt(eps0*Te/(ne0*e**2)) # Electron Debye length
LDI = np.sqrt(eps0*Ti/(ni0*e**2)) # Ion Debye Length
wpe = np.sqrt(ne0*e**2/(eps0*me)) # Electron Plasma Frequency
wpi = np.sqrt(ni0*e**2/(eps0*mi)) # Ion Plasma Frequency

mFactor = wpi/wpe
print(mFactor)

L = LDI # Characteristic Length
W = wpi # Characteristic Frequency

DT = DT_coeff*(1.0/wpe)
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data


EF = np.vstack(electric_field_data)
x = np.linspace(0,NC, EF.shape[1])
dx = x[1]-x[0]

# Combine electric field data into a 2D array
EF = EF.reshape(DATA_TS,(NC+1))
print("The shape of EF is: ", EF.shape)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wpet_1 = 0 #1000
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
E = EF[int(y1):int(y2),:]

print("The shape of E (reduced EF) is: ", E.shape)
#+++++++++++++++++++++++++ FFT of E-field data ++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define Omega and K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*L)       # Normalized k
print('Length of k: ',len(k))
print('Max of k: ',np.max(k))

Omega, K = np.meshgrid(omega, k, indexing='ij')
print('Shape of Omega: ',Omega.shape)
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
Omega /=W # Normalized Omega
K = K*L  # Normalized K
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
vnb0 = vb*(L*W) # in multiples of ion thermal velocities (v_thi)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++
raw_analytic = True
if raw_analytic:    
    wla = np.sqrt((wpe**2*(1 + 3*k**2*LDE**2)))    
    ep = wla/W    # Electron Plasma Wave (Langmuir Wave)
    kappa = k[1:]*L
    CF = 1.0/(1.0 + k[1:]**2*LDE**2) # Correction Factor due to electron
    CS = np.sqrt((Te*CF + Ti)/mi)
    w_ia = CS*k[1:] # ion acoustic mode
    #w_ia = CS*kappa # ion acoustic mode
    w_beam = vnb0*k[1:]

# The expression of the fast wave received from the paper of Gary
alpha1 = alp/(1+alp+beta)
w_fast = np.sqrt(((Ti / mi) * (1 + alpha1 + 3 * (k[1:]**2 * LDI**2) + 3 * (1 - alpha1) * (Ti / Te))) /((k[1:]* LDI)**2 + (Ti / Te) * (1 - alpha1))) * k[1:]
    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Z = np.log(np.abs(F))
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# fig, ax = plt.subplots()
# plt.pcolor(K, Omega, Z, cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
#plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=-10,vmax=5)
#plt.contourf(K, Omega, Z, cmap='rainbow')
#print(np.min(Z), np.max(Z))

extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
ax.imshow(Z, extent=extent, interpolation='nearest', aspect='auto', origin='lower')


if raw_analytic:
    ax.plot(kappa, w_ia/wpi, 'r--', linewidth = 1.0, alpha=1.0, label='IAW')
    ax.plot(kappa, w_beam/wpi, 'b--', linewidth = 1.0, alpha=1.0, label='BW')
    ax.plot(kappa, w_fast/wpi, 'm--', linewidth = 1.0, alpha=1.0, label='FW')
    

# Note: Changing vmin & vmax will change the depth of the plots.
#cbar = plt.colorbar()
#cbar.set_label('$\zeta$')

ax.set_xlabel('$k \lambda_{Di}$')
ax.set_ylabel('$\omega/\omega_{pi}$')
ax.set_xlim([0, 1.0])
ax.set_ylim([0, 3.0])
leg = ax.legend(loc='upper right',framealpha=0.5)

plt.savefig(pjoin(path,'dispersion.pdf'),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()

