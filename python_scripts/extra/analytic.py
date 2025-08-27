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
# ----------------
# --------------------------------------------------------------------------------------------
AMU = constants('atomic mass constant')
e = constants('elementary charge')
eps = constants('electric constant')

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
    # Define symbols for the coefficients
    w, k, vb ,vthe, vthi, vthn, vthb, we, wi, wn, wb = symbols('w k vb vthe vthi vthn vthb we wi wn wb')
    
    A = w**2 - k**2*vthe**2
    B = w**2 - k**2*vthi**2
    C = w**2 - k**2*vthn**2
    D = ((w - k*vb)**2) - k**2*vthb**2
    expr = A*B*C*D - we**2*B*C*D- wi**2*A*C*D - wn**2*A*B*D - wb**2*A*B*C
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
print('Highest power of omega is: %d\n' % (len(coff)-1))

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

me = constants('electron mass')

# Define the densities
ni0 = n0
#alp = 20
#beta = 0.4
ne0 = ni0 / (1 + alp + beta)
nn0 = alp*ne0
nb0 = beta * ne0
# --------------------------------------------------------------------------------------------

vthe = np.sqrt(Te / me) # Thermal velocity of electrons
vthi = np.sqrt(Ti / mi)
vthn = np.sqrt(Tn / mn)
vthb = np.sqrt(Tb / mb)

vb = vb0*vthi
print("Beam Velocity : %d vthi"%(vb))

# --------------------------------------------------------------------------------------------
we = np.sqrt((ne0 * e**2) / (eps * me))
wi = np.sqrt((ni0 * e**2) / (eps * mi))
wn = np.sqrt((nn0 * e**2) / (eps * mn))
wb = np.sqrt((nb0 * e**2) / (eps * mb))

mFactor = wi/we
print(mFactor)

LDE = np.sqrt(eps*Te/(ne0*e**2)) # Electron Debye length
LDI = np.sqrt(eps*Ti/(ni0*e**2)) # Ion Debye Length

L = LDI # Characteristic Length
W = wi # Characteristic Frequency

DT = DT_coeff*(1.0/we)

wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
#E = EF[int(y1):int(y2),:]
dx = 1
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
# -------------------------------------------------------------------------------------------
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Un-Normalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*L) 
kappa = k[1:]*LDI
CF = 1.0/(1.0 + k[1:]**2*LDE**2) # Correction Factor due to electron
CS = np.sqrt((Te*CF + Ti)/mi)
w_ia = CS*k[1:] # ion acoustic mode
#w_ia = CS*kappa # ion acoustic mode
w_beam = vb*k[1:]

# The expression of the fast wave received from the paper of Gary
alpha1 = alp/(1+alp+beta)
#beta1 = beta
#alpha1 = beta1/(1+beta1+0)
w_fast = np.sqrt(((Ti / mi) * (1 + alpha1 + 3 * (k[1:]**2 * LDI**2) + 3 * (1 - alpha1) * (Ti / Te))) /((k[1:]* LDI)**2 + (Ti / Te) * (1 - alpha1))) * k[1:]
# --------------------------------------------------------------------------------------------
# Additional parameters required
kb = constants('Boltzmann constant')
EV_TO_K = constants('electron volt-kelvin relationship')
NC = 1024
dx = 1
k = 2 * np.pi * np.arange(NC + 1) / (NC * dx * LD)

# --------------------------------------------------------------------------------------------
def EPW_dispersion():
    # Convert symbolic coefficients to numerical functions
    coeffs_func = [lambdify(['k', 'vthe', 'vthi','vthn', 'vthb','vb', 'we', 'wi','wn', 'wb'], coeff) for coeff in coff]
    
    roots = []
    for i in range(1, len(k)): 
        # Evaluate coefficients for each value of k
        coeff_values = [coeff_func(k[i], vthe, vthi,vthn,vthb,vb, we, wi,wn,wb) for coeff_func in coeffs_func]
        root = np.roots(coeff_values)
        roots.append(root)
    roots = np.array(roots)
    return roots


# --------------------------------------------------------------------------------------------
roots_EPW = EPW_dispersion()
solved_analytic = True
if solved_analytic:
    ep1 = np.real(roots_EPW[:, 0])
    ep2 = np.real(roots_EPW[:, 1])
    ep3 = np.real(roots_EPW[:, 2])
    ep4 = np.real(roots_EPW[:, 3])
    ep5 = np.real(roots_EPW[:, 4])
    ep6 = np.real(roots_EPW[:, 5])
    ep7 = np.real(roots_EPW[:, 6])
    ep8 = np.real(roots_EPW[:, 7])

    epim1 = np.imag(roots_EPW[:, 0])
    epim2 = np.imag(roots_EPW[:, 1])
    epim3 = np.imag(roots_EPW[:, 2])
    epim4 = np.imag(roots_EPW[:, 3])
    epim5 = np.imag(roots_EPW[:, 4])
    epim6 = np.imag(roots_EPW[:, 5])
    epim7 = np.imag(roots_EPW[:, 6])
    epim8 = np.imag(roots_EPW[:, 7])

    # Plot only the existing roots
    figsize = np.array([80, 80 / 1.618])  # Figure size in mm
    ppi = np.sqrt(1920**2 + 1200**2) / 24  # Screen resolution
    fig, ax = plt.subplots(1, 2, figsize=figsize / 10.4, constrained_layout=True, dpi=ppi)

    # Real part plots
    ax[0].plot(k[1:] * LD, ep1 / we, 'r.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep2 / we, 'g.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep3 / we, 'b.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep4 / we, 'k.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep5 / we, 'k.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep6 / we, 'k.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep7 / we, 'k.', markersize=1.0)
    ax[0].plot(k[1:] * LD, ep8 / we, 'k.', markersize=1.0)
    ax[0].legend()
    ax[0].set_xlabel('$k \lambda_{D}$')
    ax[0].set_ylabel('$\omega/\omega_{pe}$')
    ax[0].set_title('Real Roots')
    ax[0].grid(True)
    ax[0].set_ylim([0,5])

    # Imaginary part plots
    ax[1].plot(k[1:], epim1 / we, 'r.', markersize=2.0)
    ax[1].plot(k[1:], epim2 / we, 'g.', markersize=2.0)
    ax[1].plot(k[1:], epim3 / we, 'b.', markersize=2.0)
    ax[1].plot(k[1:], epim4 / we, 'k.', markersize=2.0)
    ax[1].plot(k[1:], epim5 / we, 'k.', markersize=2.0)
    ax[1].plot(k[1:], epim6 / we, 'k.', markersize=2.0)
    ax[1].plot(k[1:], epim7 / we, 'k.', markersize=2.0)
    ax[1].plot(k[1:], epim8 / we, 'k.', markersize=2.0)
    ax[1].set_xlabel('$k \lambda_{D}$')
    ax[1].set_ylabel('$\omega/\omega_{pe}$')
    ax[1].set_title('Imaginary Roots')
    ax[1].grid(True)

    plt.show()
