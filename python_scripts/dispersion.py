import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser
import h5py


file_name = 'result.h5'
path = sys.argv[1]

#----------------Read hdf5 file ------------
#f = h5py.File(pjoin(path, file_name), 'r', driver='core', map='read')
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']
metadata_negion = f['/metadata_species/negion']
metadata_beam = f['/metadata_species/beam']


eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

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
mi = metadata_ion.attrs['mass']
ni0 = metadata_ion.attrs['density']
nn0 = metadata_negion.attrs['density']
ne0 = metadata_electron.attrs['density']

we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']

mfactor = wp/we
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1

#DT = DT_coeff*(1.0/we)
DT = DT_coeff*(1.0/wp)

#----------------------------------------------------------------------------------
#transform elctric field data into a 2D array

electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0,NC, EF.shape[1])
dx = x[1]-x[0]
print("The shape of EF is: ", EF.shape[1])
print("The value of dx is: ", dx)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)

E = EF[int(y1):int(y2),:]
#nde = nde[int(y1):int(y2),:]

print("The shape of E (reduced EF) is: ", E.shape)
#+++++++++++++++++++++++++ FFT of E-field data ++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
#F = np.fft.fftn(nde, norm='ortho')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define Omega and K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT)*mfactor # Previously it was (NUM_TS*DT) # multiply by mfactror to convert to ion scale
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*LD) # Normalized 
print('Length of k: ',len(k))
print('Max of k: ',np.max(k))

Omega, K = np.meshgrid(omega, k, indexing='ij')
print('Shape of Omega: ',Omega.shape)
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
Omega /= wp # Normalized Omega
#Omega = Omega*mfactor
K = K*(LD)  # Normalized K
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++



Bz = 0.001

wec = e*Bz/me

wic = e*Bz/mi

omega_uh = np.sqrt(we**2 + wec**2)

omega_lh = ((wec*wic)**(-1) + wp**(-2))**(-0.5)

# Analytic EPW
raw_analytic_EPW = False
if raw_analytic_EPW:
    epe = np.sqrt(omega_uh**2 + 3 * (we * LD)**2 * k**2)

raw_analytic_upperhybrid = False
if raw_analytic_upperhybrid:
    wuh = np.full(len(k), omega_uh)

raw_analytic_lowerhybrid = False
if raw_analytic_lowerhybrid:
    wlh = np.full(len(k), omega_lh)

raw_analytic_PW = False
if raw_analytic_PW:    
    pw = np.full(len(k), we)
    #epe = epe/we


#DAW = np.sqrt(Ti*e/md)*np.sqrt(nd/1e13)
raw_analytic_IAW = False
if raw_analytic_IAW:    
    epi = np.sqrt((Te*EV_TO_K*kb/mi)*(1/(1+k**2*LD**2)) + (Ti*EV_TO_K*kb/mi))*k
    #epi = np.sqrt((Te*EV_TO_K*kb/mi)*(1/(1+k**2*LD**2))*(ni0/ne0) + (Te*EV_TO_K*kb/mi)*(nn0/ne0))*k 
    #epi = DAW*k
    #epi = wla/wp    # Electron Plasma Wave (Langmuir Wave)

alp = metadata_negion.attrs["density"]
beta = metadata_beam.attrs["density"]
alpha = alp/(1+alp+beta)
k1 = np.arange(0, 2000, 0.01)

#mi = 1e-15
raw_analytic_FIAW = True
if raw_analytic_FIAW:    
    ilw = np.sqrt(((Ti*EV_TO_K*kb/mi)*(1 + alpha + 3*(k1**2*LDi**2) + 3*(1-alpha)*(Ti/Te)))/((k1*LDi)**2 + (Ti/Te)*(1-alpha)))*k1


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#F_real = np.real(F)
#F_imag = np.imag(F)
Z = np.log(np.abs(F))
#Z = np.abs(F)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)


#-----------plotting --------------------------------------------------------

#fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
fig, ax = plt.subplots(figsize = figsize/10.4,constrained_layout=True,dpi=ppi)
c1 = plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))##
#c1 = plt.pcolor(K, Omega, Z, cmap='rainbow',shading='auto',vmin= -6,vmax = 6)
#c1 = plt.contourf(K, Omega, Z, cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
#c1 = plt.contourf(K, Omega, Z, cmap='rainbow',shading='auto',vmin=-6,vmax=6)


print(np.min(Z), np.max(Z))

#print(np.max(k*LD))

if normscheme == 2 or normscheme == 4:
    LD = LDi
    we = wp
#if normscheme == 4:
    #LD = LDi
    #we = we    
if normscheme == 1 or normscheme == 3:
    LD = LD
    we = we


print(np.max(k1*LD))


if raw_analytic_EPW:
    plt.plot(k*LD, epe/we, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$IAW$')

if raw_analytic_PW:
     plt.plot(k*LD, pw/we, color='k', linestyle='--', lw = 1.0, label='$plasma oscillation$')


if raw_analytic_IAW:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    plt.plot(k*LD, epi/we, color='k', linestyle='solid', lw = 1.0, label='$IAW$')

if raw_analytic_FIAW:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    plt.plot(k1*LD, ilw/we, color='k', linestyle='dotted', lw = 1.0, label='$FIAW$')

if raw_analytic_upperhybrid:
     plt.plot(k*LD, wuh/we, color='k', linestyle='--', lw = 1.0, label='$upper\_hybrid$')

if raw_analytic_lowerhybrid:
     plt.plot(k*LD, wlh/we, color='k', linestyle='-.', lw = 1.0, label='$lower\_hybrid$')



# Note: Changing vmin & vmax will change the depth of the plots.e
cbar = plt.colorbar()
cbar.set_label('$\zeta$')

ax.set_xlabel('$k \lambda_{De}$')
ax.set_ylabel('$\omega/\omega_{pe}$')
#ax.set_xlim([0, np.max(K)])
#ax.set_xlim([0,3])
#ax.set_ylim([0,5])
leg = ax.legend(loc='upper right',framealpha=0.5)
plt.savefig(pjoin(path,'dispersion.png'),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()

