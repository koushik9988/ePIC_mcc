import numpy as np 
from copy import deepcopy
from scipy.constants import value as constants
import matplotlib as mpp
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
from sympy import *
import sys
# --------------------------------------------------------------------------------------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
# --------------------------------------------------------------------------------------------
path = './'
# -------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
def symbolic_coeffs():
    w, k, vthi, vthe, vthn, vthnb, vnb0, wpi, wpe, wpn, wpnb = symbols('w, k, vthi, vthe, vthn, vthnb, vnb0, wpi, wpe, wpn, wpnb')
    
    A = (w**2) - (k**2)*(vthi**2)  # denominator for the background ions
    B = (w**2) - (k**2)*(vthe**2)  # denominator for the background electrons
    C = (w**2) - (k**2)*(vthn**2)  # denominator for the background negative ions 
    D = ((w - k*vnb0)**2) - (k**2)*(vthnb**2) # denominator for the beam negative ions
    
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
NUM_TS = 50000
write_interval = 5
DT_coeff = 0.01
# --------------------------------------------------------------------------------------------
Te  = 1*e
Ti  = 0.026*e
Tn  = 0.026*e
Tnb = 0.026*e
#-------------------------------------------------------
alp = 100
bet = 0.4
vb0 = 10
#------------------------------------------------------
# SIM Vars
NC = 1024
Time = 0
#------------------------------------------------------
n0 = 1E13
ni0 = n0
ne0 = n0/(1+alp+bet)
nn0 = alp*ne0
nnb0 = bet*ne0
# --------------------------------------------------------------------------------------------
# Define masses 
mi = 1*AMU  # Background ion mass
mn = 1*AMU   # Beam ion mass
mnb = 1*AMU
CS = np.sqrt(Te/mi)
#-------------------------------------------------------
DATA_TS = int(NUM_TS/write_interval) + 1
LDE = np.sqrt(eps0*Te/(ne0*e**2)) # Characteristic Debye length
LDI = np.sqrt(eps0*Ti/(ni0*e**2)) # Characteristic Debye length
# --------------------------------------------------------------------------------------------
# Define Frequencies
wpi  = np.sqrt((ni0*e**2)/(eps0*mi))    # Background ions frequency
wpe  = np.sqrt((ne0*e**2)/(eps0*me))    # Background electron Frequency
wpn = np.sqrt((nn0*e**2)/(eps0*mn))     # Beam ion frequency
wpnb = np.sqrt((nnb0*e**2)/(eps0*mnb))  # Beam electron frequency
# --------------------------------------------------------------------------------------------
DT = DT_coeff*(1.0/wpi)
# --------------------------------------------------------------------------------------------
dx = 1 # Normalized StepSize
# -------------------------------------------------------------------------------------------
# Define the polytropic coefficient
gi = 1
gnb = 1
ge = 1
gn = 1
# --------------------------------------------------------------------------------------------
# Define the Thermal Velocities
vthi = np.sqrt(gi*Ti/mi) #1E3
vthe = np.sqrt(ge*Te/me) #1E4
vthn = np.sqrt(gn*Tn/mn)
vthnb = np.sqrt(gnb*Tnb/mnb)
#--------------------------------------------------------------------------------------------
print('wpe: ', wpe)
print('wpi: ', wpi)
print('wpe/wpi:',wpe/wpi)
print('CS: ', CS)
print('vthi: ', vthi)
print('vthe: ', vthe)
#--------------------------------------------------------------------------------------------
# Define the Drift Velocities
vnb0 = vb0*vthi
#vnb0 = 90*CS 
#print('Is vthe > veb0?:',vthe>veb0)
# ---------------------------- Define k ------------------------------------------------------
wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
# -------------------------------------------------------------------------------------------
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Un-Normalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*LDI)       # un-Normalized k
# -------------------------------------------------------------------------------------------
kappa = k[1:]*LDI
CF = 1.0/(1.0 + k[1:]**2*LDE**2) # Correction Factor due to electron
CS = np.sqrt((Te*CF + gi*Ti)/mi)
w_ia = CS*k[1:] # ion acoustic mode
#w_ia = CS*kappa # ion acoustic mode
w_beam = vnb0*k[1:]

# The expression of the fast wave received from the paper of Gary
alpha1 = alp/(1+alp+bet)
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
# ---------------------------------------------------------------------------------
fig,ax = plt.subplots(1, 2, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)    
#fig, ax = plt.subplots(1, 2)
if solved_analytic:
    ms = 1.0
    ax[0].plot(kappa, (ep1/wpi), 'k-', markersize = ms)
    ax[0].plot(kappa, (ep2/wpi), 'k-', markersize = ms) 
    ax[0].plot(kappa[1:], (ep3[1:]/wpi), 'k-', markersize = ms) 
    ax[0].plot(kappa, (ep4/wpi), 'k-', markersize = ms) 
    ax[0].plot(kappa, (ep5/wpi), 'k-', markersize = ms) 
    ax[0].plot(kappa, (ep6/wpi), 'k-', markersize = ms)  
    ax[0].plot(kappa, (ep7/wpi), 'k-', markersize = ms)  
    ax[0].plot(kappa, (ep8/wpi), 'k-', markersize = ms)
    ax[0].plot(kappa[i], re[i,j], 'xr')
    #ax[0].plot(kappa[ii], ep7[ii]/wpi, 'xb')
    
    ax[0].plot(kappa, w_ia/wpi, 'r--', linewidth = 1.0, alpha=1.0, label='IAW')
    ax[0].plot(kappa, w_beam/wpi, 'b--', linewidth = 1.0, alpha=1.0, label='BW')
    ax[0].plot(kappa, w_fast/wpi, 'm--', linewidth = 1.0, alpha=1.0, label='FW')


    ax[0].set_xlabel('$k \lambda_{Di}$')
    ax[0].set_ylabel('$\omega_r/\omega_{pi}$')   
    ax[0].set_title('Real Roots')  
    ax[0].set_xlim([0, 0.5])
    ax[0].set_ylim([0, 4])
    ax[0].legend(loc=0)
    ax[0].grid(True)    
    # ------------------------------------------------------------------------------------
    print('Fastest-growing mode, k={:g}, omega={:g}'.format(kappa[i], roots_EPW[i,j]/wpi))
    # ------------------------------------------------------------------------------------
    ax[1].plot(kappa, epim1/wpi, 'k-', markersize = ms)
    ax[1].plot(kappa, epim2/wpi, 'k-', markersize = ms) 
    ax[1].plot(kappa, epim3/wpi, 'k-', markersize = ms) 
    ax[1].plot(kappa, epim4/wpi, 'k-', markersize = ms) 
    ax[1].plot(kappa[1:], epim5[1:]/wpi, 'k-', markersize = ms) 
    ax[1].plot(kappa, epim6/wpi, 'k-', markersize = ms)  
    ax[1].plot(kappa, epim7/wpi, 'k-', markersize = ms)  
    ax[1].plot(kappa, epim8/wpi, 'k-', markersize = ms)
    ax[1].plot(kappa[i], im[i,j], 'xr')
    #ax[1].plot(kappa[ii], epim7[ii]/wpi, 'xb')
    
    

    ax[1].set_xlabel('$k \lambda_{Di}$')
    ax[1].set_ylabel('$\gamma/\omega_{pi}$')   
    ax[1].set_title('Imaginary Roots')  
    ax[1].set_xlim([0, 0.5])
    #ax[1].set_ylim([0, 2])
    ax[1].grid(True)    
        

plt.savefig(pjoin(path,'dispersion_neg_beam_1.png'),dpi=dpi)   
plt.show()
