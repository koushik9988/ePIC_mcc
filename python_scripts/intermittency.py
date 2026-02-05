import numpy as np
import h5py
import sys
import os
import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.constants import value as constants
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mp


#----------------constanst--------------------

eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#----------------------------------------------

if len(sys.argv) != 3:
    print("Usage: python Dth_auto_full.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
file_name = 'result.h5'

output_dir = 'new_analysis'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

#open hdf file
f = h5py.File(pjoin(path, file_name), 'r')

#---------------------Read Metadata-----------------------------------------
metadata_group = f['/metadata']

# Read Metadata
metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
write_interval_phase = metadata.attrs['write_int_phase']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
normscheme = metadata.attrs['norm_scheme']
density = metadata.attrs['density']

species_metadata_path = f'/metadata_species/{particle_type}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']


DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

#DT = DT_coeff * (1.0 / wpi)

L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi

DT = DT_coeff * (1.0 / W)
efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12

dt_phase = write_interval_phase * DT_coeff * mfactor

# Scales to analyze
spatial_scales = np.array([1, 2, 4, 8, 16, 32, 64])
temporal_scales = np.array([1, 2, 4, 8, 16, 32])
#spatial_scales = np.array([2,6,10,14,18,22,26,30])
#temporal_scales = np.array([2,6,10,14,18,22,26,30])
n_bins = 100


electric_field_data = []

field_keys = f['fielddata/efield'].keys()
time_steps = sorted([int(k) for k in field_keys])
    
for time_step in time_steps:
    EF_data = f[f'fielddata/vel_electron/{time_step}']
    #EF_data = f[f'fielddata/den_ion/{time_step}']
    electric_field_data.append(EF_data[:])
    
EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]

print("EF shape",EF.shape)


t_start = int(100/(DT_coeff * write_interval * mfactor))#Nt // 4
t_end = Nt#3 * Nt // 4

EF_window = EF[t_start:t_end, :]
Nt_window = EF_window.shape[0]

print(f"Analyzing E(x,t): {Nt_window} timesteps x {Nx} spatial points")


spatial_flatness = []
spatial_pdfs = []


orders = np.array([1, 2, 3, 4, 5, 6])
orders = np.array([0.5,1.5, 2,2.5 ,3,3.5, 4,4.5, 5,5.5, 6])
Sp = {p: [] for p in orders} # dictionary for each p will will store arrya here [] mean fro p it currently initiazes emtyp array




for l in spatial_scales:
    # vectorized 
    """
    EF shape  = (t,x) so we take fixed time abnd vectorised space increamnet
    aand generate samples 
    """
    dE = EF_window[:, l:] - EF_window[:, :-l]

    for p in orders:
        Sp[p].append(np.mean(np.abs(dE)**p))

    dE = dE.flatten()
    dE /= np.std(dE)

    K = np.mean(dE**4)
    spatial_flatness.append(K)

    #pdf, edges = np.histogram(dE, bins=n_bins, range=(-7, 7), density=True)
    
    pdf, edges = np.histogram(dE, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    spatial_pdfs.append((centers, pdf, K))


for p in orders:
    Sp[p] = np.array(Sp[p])



temporal_flatness = []
temporal_pdfs = []


orders1 = np.array([1, 2, 3, 4, 5, 6])
orders1 = np.array([0.5,1.5, 2,2.5 ,3,3.5, 4,4.5, 5,5.5, 6])
Sp1 = {p1: [] for p1 in orders1}

for τ in temporal_scales:
    # simsilar as above
    dE = EF_window[τ:, :] - EF_window[:-τ, :]

    for p1 in orders1:
        Sp1[p1].append(np.mean(np.abs(dE)**p1))

    dE = dE.flatten()
    dE /= np.std(dE)

    K = np.mean(dE**4)
    temporal_flatness.append(K)

    pdf, edges = np.histogram(dE, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    temporal_pdfs.append((centers, pdf, K))


for p1 in orders1:
    Sp1[p1] = np.array(Sp1[p1])

dt = DT_coeff *write_interval



figsize = np.array([150,150/1.618])#Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)


fig, ax = plt.subplots(2, 2,figsize = (16,16))#,figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)

xg = np.linspace(-6, 6, 400)
Pg = np.exp(-0.5 * xg**2) / np.sqrt(2 * np.pi)

for i in range(0, len(spatial_scales), 2):
    c, p, K = spatial_pdfs[i]
    ax[0, 0].semilogy(c, p, lw=2,label=rf"$\ell={spatial_scales[i]},\ K={K:.2f}$")

ax[0, 0].semilogy(xg, Pg, 'k--', lw=2, label="Gaussian")
ax[0, 0].set_xlabel(r"$\Delta E_x / \sigma$")
ax[0, 0].set_ylabel(r"$P(\Delta E, \ell)$")
#ax.set_title("Spatial increments")
ax[0, 0].legend(fontsize=8)
#ax[0, 0].grid(True, which="both", alpha=0.3)


#ax = ax[0, 1]
for i in range(0, len(temporal_scales), 2):
    c, p, K = temporal_pdfs[i]
    ax[0,1].semilogy(c, p, lw=2,label=rf"$\tau={temporal_scales[i]},\ K={K:.2f}$")

ax[0, 1].semilogy(xg, Pg, 'k--', lw=2, label="Gaussian")
ax[0, 1].set_xlabel(r"$\Delta E_t / \sigma$")
ax[0, 1].set_ylabel(r"$P(\Delta E, \tau)$")
ax[0, 1].legend(fontsize=8)
#ax[0, 1].grid(True, which="both", alpha=0.3)


#ax = ax[1, 0]
ax[1, 0].loglog(spatial_scales * dx, spatial_flatness, 'o-', lw=2, label=r"$\ell$")
ax[1, 0].loglog(temporal_scales * dt, temporal_flatness, 's-', lw=2, label=r"$\tau$")
ax[1, 0].axhline(3, color='k', ls='--', lw=1.5, label="Gaussian")
ax[1, 0].set_xlabel("Scale")
ax[1, 0].set_ylabel("K")
ax[1, 0].legend()
#ax[1, 0].grid(True, which="both", alpha=0.3)


nt_plot = Nt_window

#ax = ax[1, 1]
im = ax[1,1].imshow(EF_window[:nt_plot],aspect='auto',origin='lower', extent=[0, NC, 
t_start*DT_coeff*write_interval*mfactor, (t_start + nt_plot)*DT_coeff*write_interval*mfactor],cmap='RdBu_r')#,vmin=-2, vmax=2)
ax[1, 1].set_xlabel("$x$")
ax[1, 1].set_ylabel(r"$\omega_{pi}t$")
#ax.set_title(r"$E(x,t)$")
plt.colorbar(im, ax=ax, label=r"$E(x,t)$")
#plt.tight_layout()

plt.savefig(pjoin(path_fig, f"intermittency.pdf"),dpi = dpi)
plt.show()


# ---- extract scaling exponents ζ_p ----
zeta = []
log_l = np.log(spatial_scales * dx)
fit_idx = slice(2, -2)

for p in orders:
    slope, _ = np.polyfit(log_l[fit_idx], np.log(Sp[p])[fit_idx], 1)
    zeta.append(slope)

zeta = np.array(zeta)

fig1, ax1 = plt.subplots(1, 3, figsize=(14, 6))

for p in orders:
    ax1[0].loglog(spatial_scales * dx, Sp[p], 'o-', lw=2, label=rf"$p={p}$")
ax1[0].set_xlabel(r"$\ell$")
ax1[0].set_ylabel(r"$S_p(\ell)$")
#ax.set_title("Structure Functions")
ax1[0].legend()
#ax1[0].grid(True, which="both", alpha=0.3)




for p1 in orders1:
    ax1[1].loglog(temporal_scales * dt, Sp1[p1], 'o-', lw=2, label=rf"$p={p1}$")
ax1[1].set_xlabel(r"$\tau$")
ax1[1].set_ylabel(r"$S_p(\tau)$")
#ax.set_title("Structure Functions")
ax1[1].legend()


ax1[2].plot(orders, zeta, 'o-', lw=2, label=r"$\zeta_p$")
ax1[2].plot(orders, orders/3, 'k--', lw=2, label=r"K41")
ax1[2].set_xlabel(r"$p$")
ax1[2].set_ylabel(r"$\zeta_p$")
#ax.set_title("Scaling Exponents")
ax1[2].legend()
#ax1[1].grid(True, alpha=0.3)

#ax = axs2[2]
#x = np.asarray(x)
#ax.plot(x**(2/3),y)


plt.tight_layout()
plt.savefig(pjoin(path_fig, "structure_function_scaling.pdf"), dpi=dpi)
plt.show()
