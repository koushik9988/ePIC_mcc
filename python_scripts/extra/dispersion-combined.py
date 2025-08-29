import numpy as np 
from copy import deepcopy
from scipy.constants import value as constants
import matplotlib as mpp
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
from sympy import *
import sys
import h5py
# --------------------------------------------------------------------------------------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
# --------------------------------------------------------------------------------------------
file_name = 'result.h5'
#path = '../data_run_3/'
path = sys.argv[1]
f = h5py.File(pjoin(path, file_name), 'r')
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']
metadata_negion = f['/metadata_species/negion']
metadata_beam = f['/metadata_species/beam']
metadata_group = f['/metadata']

# -------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
def symbolic_coeffs():
    w, k, vthi, vthe, vthn, vthb, vnb0, wpi, wpe, wpn, wpnb = symbols('w, k, vthi, vthe, vthn, vthb, vnb0, wpi, wpe, wpn, wpnb')
    
    A = (w**2) - (k**2)*(vthi**2)  # denominator for the background ions
    B = (w**2) - (k**2)*(vthe**2)  # denominator for the background electrons
    C = (w**2) - (k**2)*(vthn**2)  # denominator for the background negative ions 
    D = ((w - k*vnb0)**2) - (k**2)*(vthb**2) # denominator for the beam negative ions
    
    expr =  A*B*C*D - (B*C*D*wpi**2) - (A*C*D*wpe**2) - (A*B*D*wpn**2) - (A*B*C*wpnb**2)
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
# np.roots doesn't always arrange the roots in the same order.
# This function is used to correct that.
# Consecutive samples are reordered to have the smallest possible jump.
# Function written by: Dr. Sigvald Marholm
def correct(roots):
    nSamples, nRoots = roots.shape
    for i in range(1,nSamples):
        for j in range(nRoots):
            dist = np.abs(roots[i,:]-roots[i-1,j])
            closest = np.argmin(dist)
            if closest != j:
                # Swap root j and closest
                tmp = deepcopy(roots[i:,j])
                roots[i:,j] = roots[i:,closest]
                roots[i:,closest] = tmp
                #print('SWAP')
    return roots
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
#print(coff)
print('Highest power of omega is: %d\n'%(len(coff)-1))
# -------------------------------------------------------------------------------------------
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
vb0 = metadata_beam.attrs['streaming_velocity']
nParticlesE = metadata_electron.attrs['num_particles']
nParticlesI = metadata_ion.attrs['num_particles']
nParticlesN = metadata_negion.attrs['num_particles']
nParticlesB = metadata_beam.attrs['num_particles']
alp = metadata_negion.attrs['density']
beta = metadata_beam.attrs['density']

EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1
#------------------------------------------------------
vthi = np.sqrt(Ti/mi)
vthe = np.sqrt(Te/me)
vthn = np.sqrt(Tn/mi)
vthb = np.sqrt(Tb/mi)
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
wpn = np.sqrt((nn0*e**2)/(eps0*mn))     # Beam ion frequency
wpnb = np.sqrt((nb0*e**2)/(eps0*mb))  # Beam electron frequency

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
EF = EF.reshape(DATA_TS,(NC+0))
print("The shape of EF is: ", EF.shape)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the Drift Velocities
vnb0 = vb0*(L*W)
print("Beam Velocity : %d vthi"%(vb0))
# ---------------------------- Define k ------------------------------------------------------
wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
E = EF[int(y1):int(y2),:]
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
# -------------------------------------------------------------------------------------------
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Un-Normalized Omega
k     = 2*np.pi*np.arange(NC+0)/(NC*dx*L)       # un-Normalized k
Omega, K = np.meshgrid(omega, k, indexing='ij')
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
Omega /=W # Normalized Omega
K = K*L  # Normalized K
# -------------------------------------------------------------------------------------------
kappa = k[1:]*LDI
CF = 1.0/(1.0 + k[1:]**2*LDE**2) # Correction Factor due to electron
CS = np.sqrt((Te*CF + Ti)/mi)
w_ia = CS*k[1:] # ion acoustic mode
#w_ia = CS*kappa # ion acoustic mode
w_beam = vnb0*k[1:]

# The expression of the fast wave received from the paper of Gary
alpha1 = alp/(1+alp+beta)
#beta1 = beta
#alpha1 = beta1/(1+beta1+0)
w_fast = np.sqrt(((Ti / mi) * (1 + alpha1 + 3 * (k[1:]**2 * LDI**2) + 3 * (1 - alpha1) * (Ti / Te))) /((k[1:]* LDI)**2 + (Ti / Te) * (1 - alpha1))) * k[1:]
# --------------------------------------------------------------------------------------------
#                       ANALYTICAL CALCULATIONS
# --------------------------------------------------------------------------------------------

def EPW_dispersion():
    # Coefficients
     coeff1 = eval(str(coff[0]))
     coeff2 = eval(str(coff[1]))
     coeff3 = eval(str(coff[2]))
     coeff4 = eval(str(coff[3]))
     coeff5 = eval(str(coff[4]))
     coeff6 = eval(str(coff[5]))
     coeff7 = eval(str(coff[6]))
     coeff8 = eval(str(coff[7]))
     coeff9 = eval(str(coff[8]))
          
     roots = []
     for i in range(1,len(k)): 
         coeffs = [coeff1, coeff2[i], coeff3[i], coeff4[i], coeff5[i], coeff6[i], coeff7[i], coeff8[i], coeff9[i]]
         root = np.roots(coeffs)
         roots.append(root)
     roots = np.array(roots)
     roots = correct(roots)
     return roots
# --------------------------------------------------------------------------------------------
roots_EPW = EPW_dispersion()
solved_analytic = True
if solved_analytic:    
    ep1 = np.real(roots_EPW[:,0])
    ep2 = np.real(roots_EPW[:,1])
    ep3 = np.real(roots_EPW[:,2])
    ep4 = np.real(roots_EPW[:,3])
    ep5 = np.real(roots_EPW[:,4])
    ep6 = np.real(roots_EPW[:,5])
    ep7 = np.real(roots_EPW[:,6])
    ep8 = np.real(roots_EPW[:,7])    
    #----------------------------------
    epim1 = np.imag(roots_EPW[:,0])
    epim2 = np.imag(roots_EPW[:,1])
    epim3 = np.imag(roots_EPW[:,2])
    epim4 = np.imag(roots_EPW[:,3])
    epim5 = np.imag(roots_EPW[:,4])
    epim6 = np.imag(roots_EPW[:,5])
    epim7 = np.imag(roots_EPW[:,6])
    epim8 = np.imag(roots_EPW[:,7])
        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([100,100/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)
# ---------------------------------------------------------------------------------
re = np.real(roots_EPW/wpi) # all real parts of the roots
im = np.imag(roots_EPW/wpi) # all imaginary parts of the roots
index = np.argmax(im)
print('index', index)
i,j = np.unravel_index(index, im.shape)
# ---------------------------------------------------------------------------------
ep7im =  epim7/wpi
ind = np.argmax(ep7im)
ii = np.unravel_index(ind, ep7im.shape)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Z = np.log(np.abs(F))
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#               FOR BOTH REAL & IMAGINARY PLOT TOGETHER
# ---------------------------------------------------------------------------- 
# fig,ax = plt.subplots(1, 2, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)    

# if solved_analytic:
#     extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
#     ax[0].imshow(Z, extent=extent, interpolation='nearest', aspect='auto', origin='lower')

#     ms = 0.8
#     ax[0].plot(kappa, (ep1/wpi), 'k-', markersize = ms, linewidth = ms)
#     ax[0].plot(kappa, (ep2/wpi), 'k-', markersize = ms, linewidth = ms) 
#     ax[0].plot(kappa[1:], (ep3[1:]/wpi), 'k-', markersize = ms, linewidth = ms) 
#     ax[0].plot(kappa, (ep4/wpi), 'k-', markersize = ms, linewidth = ms) 
#     ax[0].plot(kappa, (ep5/wpi), 'k-', markersize = ms, linewidth = ms) 
#     ax[0].plot(kappa, (ep6/wpi), 'k-', markersize = ms, linewidth = ms)  
#     ax[0].plot(kappa, (ep7/wpi), 'k-', markersize = ms, linewidth = ms)  
#     ax[0].plot(kappa, (ep8/wpi), 'k-', markersize = ms, linewidth = ms)
#     ax[0].plot(kappa[i], re[i,j], 'xr')
#     #ax[0].plot(kappa[ii], ep7[ii]/wpi, 'xb')
    
#     ax[0].plot(kappa, w_ia/wpi, 'r--', linewidth = 0.5, alpha=0.5, label='IAW')
#     ax[0].plot(kappa, w_beam/wpi, 'b--', linewidth = 0.5, alpha=0.5, label='BW')
#     ax[0].plot(kappa, w_fast/wpi, 'm--', linewidth = 0.5, alpha=0.5, label='FW')


#     ax[0].set_xlabel('$k \lambda_{Di}$')
#     ax[0].set_ylabel('$\omega_r/\omega_{pi}$')   
#     ax[0].set_title('Real Roots')  
#     ax[0].set_xlim([0, 1.0])
#     ax[0].set_ylim([0, 4])
#     ax[0].legend(loc=0, framealpha=0.5)
#     ax[0].grid(True)    
#     # ------------------------------------------------------------------------------------
#     print('Fastest-growing mode, k={:g}, omega={:g}'.format(kappa[i], roots_EPW[i,j]/wpi))
#     # ------------------------------------------------------------------------------------
#     ax[1].plot(kappa, epim1/wpi, 'k-', markersize = ms)
#     ax[1].plot(kappa, epim2/wpi, 'k-', markersize = ms) 
#     ax[1].plot(kappa, epim3/wpi, 'k-', markersize = ms) 
#     ax[1].plot(kappa, epim4/wpi, 'k-', markersize = ms) 
#     ax[1].plot(kappa[1:], epim5[1:]/wpi, 'k-', markersize = ms) 
#     ax[1].plot(kappa, epim6/wpi, 'k-', markersize = ms)  
#     ax[1].plot(kappa, epim7/wpi, 'k-', markersize = ms)  
#     ax[1].plot(kappa, epim8/wpi, 'k-', markersize = ms)
#     ax[1].plot(kappa[i], im[i,j], 'xr')
#     #ax[1].plot(kappa[ii], epim7[ii]/wpi, 'xb')  
    

#     ax[1].set_xlabel('$k \lambda_{Di}$')
#     ax[1].set_ylabel('$\gamma/\omega_{pi}$')   
#     ax[1].set_title('Imaginary Roots')  
#     ax[1].set_xlim([0, 1.0])
#     #ax[1].set_ylim([0, 2])
#     ax[1].grid(True)    
        
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                       FOR REAL PLOT & HIGHLIGHT
# ---------------------------------------------------------------------------- 
fig,ax = plt.subplots(1, 1, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)    

if solved_analytic:
    extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
    ax.imshow(Z, extent=extent, interpolation='nearest', aspect='auto', origin='lower')

    ms = 0.8
    ax.plot(kappa, (ep1/wpi), 'k-', markersize = ms, linewidth = ms)
    ax.plot(kappa, (ep2/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax.plot(kappa[1:], (ep3[1:]/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax.plot(kappa, (ep4/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax.plot(kappa, (ep5/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax.plot(kappa, (ep6/wpi), 'k-', markersize = ms, linewidth = ms)  
    ax.plot(kappa, (ep7/wpi), 'k-', markersize = ms, linewidth = ms)  
    ax.plot(kappa, (ep8/wpi), 'k-', markersize = ms, linewidth = ms)
    ax.plot(kappa[i], re[i,j], 'xr')
    #ax[0].plot(kappa[ii], ep7[ii]/wpi, 'xb')
    
    ax.plot(kappa, w_ia/wpi, 'r--', linewidth = 0.5, alpha=0.5, label='IAW')
    ax.plot(kappa, w_beam/wpi, 'b--', linewidth = 0.5, alpha=0.5, label='BW')
    ax.plot(kappa, w_fast/wpi, 'm--', linewidth = 0.5, alpha=0.5, label='FW')

     
    # ------------------------------------------------------------------------------------
    print('Fastest-growing mode, k={:g}, omega={:g}'.format(kappa[i], roots_EPW[i,j]/wpi))
    # ------------------------------------------------------------------------------------
    ax.set_xlabel('$k \lambda_{Di}$')
    ax.set_ylabel('$\omega_r/\omega_{pi}$')   
    # ax.set_title('Real Roots')  
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0, 4])
    ax.legend(loc=0, framealpha=0.5)
    ax.grid(False)

    # ------- inset plot -------------------
    # set the coordinates of the main plot to be zoomed
    x1, x2, y1, y2 = 0.02, 0.18, 0.7, 1.8 
    ax_inset = ax.inset_axes([0.3, 0.5, 0.3, 0.3],
                        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    ax_inset.imshow(Z, extent=extent, origin="lower", interpolation='bilinear', aspect='auto')
    ax_inset.plot(kappa, w_fast/wpi, 'm--', linewidth = 0.5, alpha=0.5, label='FW')
    ax_inset.plot(kappa, (ep1/wpi), 'k-', markersize = ms, linewidth = ms)
    ax_inset.plot(kappa, (ep2/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax_inset.plot(kappa[1:], (ep3[1:]/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax_inset.plot(kappa, (ep4/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax_inset.plot(kappa, (ep5/wpi), 'k-', markersize = ms, linewidth = ms) 
    ax_inset.plot(kappa, (ep6/wpi), 'k-', markersize = ms, linewidth = ms)  
    ax_inset.plot(kappa, (ep7/wpi), 'k-', markersize = ms, linewidth = ms)  
    ax_inset.plot(kappa, (ep8/wpi), 'k-', markersize = ms, linewidth = ms)
    ax_inset.plot(kappa[i], re[i,j], 'xr')
    ax_inset.plot(kappa, w_beam/wpi, 'b--', linewidth = ms, alpha = 0.5)
    
    ax.plot(kappa, w_ia/wpi, 'r--', linewidth = 0.5, alpha=0.5, label='IAW')
    ax.plot(kappa, w_beam/wpi, 'b--', linewidth = 0.5, alpha=0.5, label='BW')
    ax.plot(kappa, w_fast/wpi, 'm--', linewidth = 0.5, alpha=0.5, label='FW')
    
    ax.indicate_inset_zoom(ax_inset, facecolor = "none", edgecolor="blue", linewidth=1.5, alpha=0.5, zorder=4.99)

   


plt.savefig(pjoin(path,'dispersion_combined.pdf'),dpi=dpi)   
plt.show()
