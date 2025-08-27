import numpy as np 
from scipy.constants import value as constants
import matplotlib as mpp
import matplotlib.pyplot as plt
from sympy import *
# --------------------------------------------------------------------------------------------
AMU = constants('atomic mass constant')
e = constants('elementary charge')
eps = constants('electric constant')
# --------------------------------------------------------------------------------------------
def symbolic_coeffs():
    #w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn = symbols('w, k, v10, v20, ve0, vth1, vth2, vthe, w1, w2, we, wn')
    w, k, vb ,vthe, vthi, vthn, vthb, we, wi, wn, wb = symbols('w, k, vb ,vthe, vthi, vthn, vthb, we, wi, wn, wb')
    
    A = w**2 - k**2*vthe**2
    B = w**2 - k**2*vthi**2
    C = w**2 - k**2*vthn**2
    D = ((w - k*vb)**2) - k**2*vthb**2
    expr =  A*B*C*D - we**2*B*C*D - wi**2*A*C*D - wn**2*A*B*D - wb**2*A*B*C
    #expr = simplify(expr)  # Add this line after defining expr
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
#print(coff)
print('Highest power of omega is: %d\n'%(len(coff)-1))
# --------------------------------------------------------------------------------------------
# Define the densities
ni0 = 1e13
alpha = 20
beta = 0.4
ne0 = ni0/(1+alpha+beta)
nn0 = alpha*ne0
nb0 = beta*ne0
#---------------------------------------------------------------------------------------------
# Define masses 
me = constants('electron mass')
mi = 1*AMU
mb = 1*AMU
mn = 1*AMU 

Te = 1*e
Ti = 0.026*e

vthe = np.sqrt(Te/me) #1E4
vthi = np.sqrt(Ti/mi)
vthb = vthi
vthn = vthi

vb = 10*vthi

k = np.arange(0, 500, 0.01)
# --------------------------------------------------------------------------------------------
we = np.sqrt((ne0*e**2)/(eps*me))
wn = np.sqrt((nn0*e**2)/(eps*mn))
wi = np.sqrt((ni0*e**2)/(eps*mn))
wb = np.sqrt((nb0*e**2)/(eps*mn))
# --------------------------------------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution



#alpha1 for analytic formula for fast mode
alpha1 = alpha / (1 + alpha + beta)
#---------------------------------------------------------------------------------------------
# Additional parameters required for ILW
kb = constants('Boltzmann constant')
EV_TO_K = constants('electron volt-kelvin relationship')
LDi = vthi / wi  # Debye length for ions, assuming normalization with electron frequency we
LD = vthe/we

#LDi = np.sqrt(eps*Ti/ni0*e*e)
NC = 1024
dx = 1 
k   = 2*np.pi*np.arange(NC+1)/(NC*dx*LD) 
#---------------------------------------------------------------------------------------------

iaw = np.sqrt(((Ti / mi) * (1 + alpha1 + 3 * (k**2 * LDi**2) + 3 * (1 - alpha1) * (Ti / Te))) /((k * LDi)**2 + (Ti / Te) * (1 - alpha1))) * k

iaw1 = np.sqrt((Te/mi)*(1/(1+k**2*LD**2)) + (Ti*EV_TO_K*kb/mi))*k

k_plot = k[1:]  
iaw_plot = iaw[1:]  # Slice ilw to match k_plot
iaw1_plot = iaw1[1:]  # Slice ilw to match k_plot

epe = np.sqrt(we**2 + 2*(3/2)*(we*LD)**2*k**2)

epwplot = epe[1:]
#############################
#LD = LDi
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
     print(coeff2)

     roots = []
     for i in range(1,len(k)): 
         coeffs = [coeff1, coeff2[i], coeff3[i], coeff4[i], coeff5[i], coeff6[i], coeff7[i], coeff8[i],coeff9[i]]
         root = np.roots(coeffs)
         roots.append(root)
     roots = np.array(roots)
     return roots
# -------------------


roots_EPW = EPW_dispersion()

#print(roots_EPW[:,1])
#print(roots_EPW[:,3])
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)

solved_analytic = True
LD = LDi
we = wi
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if solved_analytic:

    ep1 = np.real(roots_EPW[:,0])
    ep2 = np.real(roots_EPW[:,1])
    ep3 = np.real(roots_EPW[:,2])
    ep4 = np.real(roots_EPW[:,3])
    ep5 = np.real(roots_EPW[:,4])
    ep6 = np.real(roots_EPW[:,5])
    ep7 = np.real(roots_EPW[:,6])
    ep8 = np.real(roots_EPW[:,7])

    epim1 = np.imag(roots_EPW[:,0])
    epim2 = np.imag(roots_EPW[:,1])
    epim3 = np.imag(roots_EPW[:,2])
    epim4 = np.imag(roots_EPW[:,3])
    epim5 = np.imag(roots_EPW[:,4])
    epim6 = np.imag(roots_EPW[:,5])
    epim7 = np.imag(roots_EPW[:,6])
    epim8 = np.imag(roots_EPW[:,7])

    fig, ax = plt.subplots(1, 2, figsize=figsize/10.4,constrained_layout=True,dpi=ppi)
    #ax[0].semilogx(k[1:], ep1/we, color='r', linestyle='-', lw = 2.0, label='$EPW$')
    #ax[0].semilogx(k[1:], ep2/we, color='g', linestyle='-', lw = 2.0, label='$EPW$') 
    #ax[0].semilogx(k[1:], ep3/we, color='b', linestyle='-', lw = 2.0, label='$EPW$') 
    #ax[0].semilogx(k[1:], ep4/we, color='k', linestyle='-', lw = 2.0, label='$EPW$') 
    #ax[0].semilogx(k[1:], ep5/we, color='c', linestyle='-', lw = 2.0, label='$EPW$') 
    #ax[0].semilogx(k[1:], ep6/we, color='m', linestyle='-', lw = 2.0, label='$EPW$')

    ax[0].plot(k[1:]*LD, ep1/we, 'r.', markersize = 1.0)
    ax[0].plot(k[1:]*LD, ep2/we, 'g.', markersize = 1.0) 
    ax[0].plot(k[1:]*LD, ep3/we, 'b.', markersize = 1.0) 
    ax[0].plot(k[1:]*LD, ep4/we, 'k.', markersize = 1.0) 
    ax[0].plot(k[1:]*LD, ep5/we, 'c.', markersize = 1.0) 
    ax[0].plot(k[1:]*LD, ep6/we, 'm.', markersize = 1.0)  
    ax[0].plot(k[1:]*LD, ep7/we, 'y.', markersize = 1.0)  
    ax[0].plot(k[1:]*LD, ep8/we, 'orange', markersize = 1.0) 
    ax[0].plot(k[1:]*LD, iaw_plot/we, color='k', linestyle='-.', lw=1.0, label='$FIAW$')
    ax[0].plot(k[1:]*LD, iaw1_plot/we, color='k', linestyle='-.', lw=1.0, label='$IAW$') 
    
    ax[0].set_xlabel('$k \lambda_{D}$')
    ax[0].set_ylabel('$\omega/\omega_{pe}$')   
    ax[0].set_title('Real Roots')
    ax[0].grid(True)
    #ax[0].set_xlim([0, np.max(k)])
    ax[0].set_ylim([0, 2])
    ax[0].legend()

    ax[1].plot(k[1:]*LD, epim1/we, 'r.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim2/we, 'g.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim3/we, 'b.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim4/we, 'k.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim5/we, 'c.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim6/we, 'm.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim7/we, 'y.', markersize = 2.0)
    ax[1].plot(k[1:]*LD, epim8/we, 'orange', markersize = 2.0)
    
    ax[1].set_xlabel('$k \lambda_{D}$')
    ax[1].set_ylabel('$\omega/\omega_{pe}$')  
    ax[1].set_title('Imaginary Roots') 
    ax[1].grid(True)

plt.show()
