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
    # Define symbols for the coefficients
    w, k, vb ,vthe, vthi, vthb, we, wi, wb = symbols('w k vb vthe vthi vthb we wi wb')
    
    A = w**2 - k**2*vthe**2
    B = w**2 - k**2*vthi**2
    C = ((w - k*vb)**2) - k**2*vthb**2
    expr = A*B*C - we**2*B*C - wi**2*A*C - wb**2*A*B
    p = Poly(expr, w)
    coff = p.all_coeffs()    
    return coff
# --------------------------------------------------------------------------------------------
coff = symbolic_coeffs()
print('Highest power of omega is: %d\n' % (len(coff)-1))

# Define the densities
ni0 = 1e13
beta = 10
ne0 = ni0 / (1 + beta)
nb0 = beta * ne0
# --------------------------------------------------------------------------------------------
# Define masses
me = constants('electron mass')
mi = 1 * AMU
mb = constants('electron mass')

Te = 1 * e
Ti = 3* e
Tb = 0.1* e

vthe = np.sqrt(Te / me) # Thermal velocity of electrons
vthi = np.sqrt(Ti / mi)
vthb = np.sqrt(Tb / me)

vb = 5*vthe

# --------------------------------------------------------------------------------------------
we = np.sqrt((ne0 * e**2) / (eps * me))
wi = np.sqrt((ni0 * e**2) / (eps * mi))
wb = np.sqrt((nb0 * e**2) / (eps * me))

#LDi = vthi / wi  # Debye length for ions
LD = vthe / we

# --------------------------------------------------------------------------------------------
# Additional parameters required
kb = constants('Boltzmann constant')
EV_TO_K = constants('electron volt-kelvin relationship')
NC = 1024
dx = 1
k = 2 * np.pi * np.arange(NC + 1) / (NC * dx * LD)
k = np.linspace(0,5000,10000)

# --------------------------------------------------------------------------------------------
def EPW_dispersion():
    # Convert symbolic coefficients to numerical functions
    coeffs_func = [lambdify(['k', 'vthe', 'vthi','vthb', 'vb', 'we', 'wi', 'wb'], coeff) for coeff in coff]
    
    roots = []
    for i in range(1, len(k)): 
        # Evaluate coefficients for each value of k
        coeff_values = [coeff_func(k[i], vthe, vthi,vthb,vb, we, wi, wb) for coeff_func in coeffs_func]
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

    epim1 = np.imag(roots_EPW[:, 0])
    epim2 = np.imag(roots_EPW[:, 1])
    epim3 = np.imag(roots_EPW[:, 2])
    epim4 = np.imag(roots_EPW[:, 3])
    epim5 = np.imag(roots_EPW[:, 4])
    epim6 = np.imag(roots_EPW[:, 5])

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
    ax[1].set_xlabel('$k \lambda_{D}$')
    ax[1].set_ylabel('$\omega/\omega_{pe}$')
    ax[1].set_title('Imaginary Roots')
    ax[1].grid(True)

    plt.show()
