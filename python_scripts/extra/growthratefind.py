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
#print('Highest power of omega is: %d\n'%(len(coff)-1))
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
vb = metadata_beam.attrs['streaming_velocity']
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
#print(mFactor)

L = LDI # Characteristic Length
W = wpi # Characteristic Frequency

DT = DT_coeff*(1.0/wpe)

#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++

# Define the Drift Velocities
vnb0 = vb*(L*W)
print("Beam Velocity : %d vthi"%(vb))
dx = 1
wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
# -------------------------------------------------------------------------------------------
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Un-Normalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*L)       # un-Normalized k

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
dpi = 1200                        #Print resolution
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
#print('index', index)
i,j = np.unravel_index(index, im.shape)
# ---------------------------------------------------------------------------------
ep7im =  epim7/wpi
ind = np.argmax(ep7im)
ii = np.unravel_index(ind, ep7im.shape)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Transform electric field data into a 2D array
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for p, time_step in enumerate(time_steps):
    EF_data = f['fielddata/efield/' + str(time_step)]
    electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0, NC, EF.shape[1], endpoint=False)  # Spatial grid
dx = x[1] - x[0]
#print("The shape of EF is: ", EF.shape)
#print("The value of dx is: ", dx)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select time range 
wpet_1 = 0 
wpet_2 = NUM_TS * DT_coeff *0.4
y1 = int(wpet_1 / (DT_coeff * write_interval))
y2 = int(wpet_2 / (DT_coeff * write_interval))

E = EF[y1:y2, :]  # Electric field in linear phase
#print("Shape of E (linear phase): ", E.shape)

#+++++++++++++++++++++++++ Spatial FFT to Find k ++++++++++++++++++++++++++++++++
# Perform FFT along spatial axis (axis=1) for each timestep (1d fft done along spatial axis over time (if axis = 0 then 1d fft along time axis done for x space))
F_k = np.fft.fft(E, axis=1, norm='ortho')  # Shape: (time_steps, k_modes)
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)  # Wavenumber array (unnormalized)
#k = k * LDi # Normalize k by Debye length

#power spectrum
power_spectrum = np.mean(np.abs(F_k)**2, axis=0)

# Find dominant k
k_max_idx = np.argmax(power_spectrum[:NC//2])
k_max = kfft[k_max_idx]
#print(f"Dominant k (normalized by LD): {k_max:.4f}")


# Plot power spectrum
plt.figure()
plt.plot(kfft[:NC//2], power_spectrum[:NC//2])
plt.axvline(k_max, color='r', linestyle='--', label=f'k_max = {k_max:.4f}')
plt.xlabel('k')
plt.ylabel('$|E(k)|^2$')
plt.title('Power Spectrum')
#plt.legend()
plt.grid()
plt.show()

time = np.arange(y1, y2)*DT*write_interval  # Time array unnormalized 
E_k_max = np.abs(F_k[:, k_max_idx])  # Amplitude of dominant k over time

#fit_range = slice(0, len(time)//2)  # Adjust to fit linear phase
fit_range = slice(0, len(time))  # Adjust to fit linear phase

coeffs = np.polyfit(time[fit_range], np.log(E_k_max[fit_range]), 1)
gamma = coeffs[0]/wpi
#plt.semilogy(time, np.exp(coeffs[1] + coeffs[0]*time), 'r--', label=f'Fit: γ = {gamma:.2e}')
#plt.legend()
#plt.show()


#fig,(ax,ax1) = plt.subplots(1, 2, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)    
fig,(ax,ax1) = plt.subplots(1, 2) 

if solved_analytic:

    ax.semilogy(time, E_k_max)
    ax.semilogy(time, np.exp(coeffs[1] + coeffs[0]*time), 'r--', label=f'$\gamma/\omega_{{pi}} = {gamma:.4f}$')
    ax.set_xlabel("time")
    ax.set_ylabel('$\log E$')
    #ax.set_title("Plot of $\log E$ vs time to Find growth rate")
    ax.legend()

    ms = 0.8
    ax1.plot(kappa, epim1/wpi, 'k-', markersize = ms)
    ax1.plot(kappa, epim2/wpi, 'k-', markersize = ms) 
    ax1.plot(kappa, epim3/wpi, 'k-', markersize = ms) 
    ax1.plot(kappa, epim4/wpi, 'k-', markersize = ms) 
    ax1.plot(kappa[1:], epim5[1:]/wpi, 'k-', markersize = ms) 
    ax1.plot(kappa, epim6/wpi, 'k-', markersize = ms)  
    ax1.plot(kappa, epim7/wpi, 'k-', markersize = ms)  
    ax1.plot(kappa, epim8/wpi, 'k-', markersize = ms)

    #ax1.plot(kappa[i], im[i,j], 'xr')
    #ax1.plot(kappa[ii], epim7[ii]/wpi, 'xb')


    k_fast = kappa[i]
    #k_fast = k_max
    gamma_fast = np.imag(roots_EPW[i,j]/wpi)

    ax1.plot(kappa[i], im[i,j], 'xr',label=f'k = {kappa[i]:.2f}, $\gamma/\omega_{{pi}} = {np.imag(roots_EPW[i,j]/wpi):.4f}$')
    #ax1.plot(kappa[ii], epim7[ii]/wpi, 'xb')  

    # Draw lines from axes to the marker
    ax1.axvline(x=k_fast, ymin=0, ymax=gamma_fast/0.20, color='r', ls='--', lw=0.5, alpha=0.5)  # Vertical from x-axis
    ax1.axhline(y=gamma_fast, xmin=0, xmax=k_fast/1, color='r', ls='--', lw=0.5, alpha=0.5)  # Horizontal from y-axis
    

    ax1.set_xlabel('$k \lambda_{Di}$')
    ax1.set_ylabel('$\gamma/\omega_{pi}$')   
    ##ax1.set_title('Plot of imaginary part of analytical disperison relation')  
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0,0.20])
    ax1.grid(True) 
    ax1.legend()

    # ------------------------------------------------------------------------------------
    print('Fastest-growing mode, k={:g}, omega={:g}'.format(kappa[i], roots_EPW[i,j]/wpi))
    print(f"Growth rate γ for k from PIC = {k_max:.4f}: {gamma:.4f}")

plt.savefig(pjoin(path,'growth_test.pdf'),dpi=dpi)   
plt.show()
