import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from os.path import join as pjoin
import matplotlib as mp

file_name = 'result.h5'

# Fixed paths (replace with your sys.argv if needed)
path1 = "../datafiles/data_detach_1e18_vb_10"
path2 = "../datafiles/data_detach_1e18_vb_15"
path3 = "../datafiles/data_detach_1e18_vb_20"
path4 = "../datafiles/data_detach_1e18"
path5 = "../datafiles/data_detach_1e18_vb_30"

output_dir = './plots'
path_fig = pjoin(path1, output_dir)
os.makedirs(path_fig, exist_ok=True)

# Open HDF5 from first path to get time
f = h5py.File(pjoin(path1, file_name), 'r')
metadata_group = f['/metadata']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
mfactor = wpi / wpe
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor  # Time
f.close()

def load_data(path):
    f = h5py.File(pjoin(path, file_name), 'r')
    vd = f['/metadata_species/beam'].attrs['streaming_velocity']
    data = f["time_var/kinetic_energy"]
    data1 = f["time_var/electronegativity"]
    spno = f["/metadata"].attrs["spno"]
    pe = data[:, 3 * spno + 1]  # Potential energy column
    alpha = data1[:, 1]  # Electronegativity column
    f.close()
    return pe, alpha, vd


def density_to_pressure(density, T_ev=0.026):
    kB = 1.380649e-23  # J/K
    eV_to_K = 11604.525
    T = T_ev * eV_to_K  # in K
    P_Pa = density * kB * T
    P_torr = P_Pa / 133.322
    return P_torr

p = density_to_pressure(1e18)
print(p)

# Load data for all 5 paths
pe1, alpha1, vd1 = load_data(path1)
pe2, alpha2, vd2 = load_data(path2)
pe3, alpha3, vd3 = load_data(path3)
pe4, alpha4, vd4 = load_data(path4)
pe5, alpha5, vd5 = load_data(path5)

# ---------- Plot 1: Time vs Potential Energy ----------
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=15, serif='Computer Modern Roman')
mp.rc('axes', titlesize=15)
mp.rc('axes', labelsize=15)
mp.rc('xtick', labelsize=15)
mp.rc('ytick', labelsize=15)
mp.rc('legend', fontsize=15)


fig1, ax1 = plt.subplots(figsize=figsize/10.4, constrained_layout=True, dpi=ppi)

ax1.plot(ts, pe1, label=fr"v$_d$ = {int(vd1)}", color='purple')
ax1.plot(ts, pe2, label=fr"v$_d$ = {int(vd2)}", color='blue')
ax1.plot(ts, pe3, label=fr"v$_d$ = {int(vd3)}", color='green')
ax1.plot(ts, pe4, label=fr"v$_d$ = {int(vd4)}", color='red')
ax1.plot(ts, pe5, label=fr"v$_d$ = {int(vd5)}", color='orange')
ax1.set_ylim([1e-4,1e1])

ax1.set_xlabel(r"$\omega_{pi} t$")
ax1.set_ylabel("Potential Energy")
ax1.semilogy()
ax1.set_title(f"Time Evolution of Potential Energy at Pressure = {p:.2e} torr")
#ax1.grid(True, which='both', linestyle='--', alpha=0.5)
ax1.legend()
plt.tight_layout()
plt.savefig(pjoin(path_fig, 'time_vs_pe.pdf'), dpi=1200)

# ---------- Plot 2: Time vs Electronegativity ----------
fig2, ax2 = plt.subplots(figsize=figsize/10.4, constrained_layout=True, dpi=ppi)

ax2.plot(ts, alpha1, label=fr"v$_d$ = {int(vd1)}", color='purple')
ax2.plot(ts, alpha2, label=fr"v$_d$ = {int(vd2)}", color='blue')
ax2.plot(ts, alpha3, label=fr"v$_d$ = {int(vd3)}", color='green')
ax2.plot(ts, alpha4, label=fr"v$_d$ = {int(vd4)}", color='red')
ax2.plot(ts, alpha5, label=fr"v$_d$ = {int(vd5)}", color='orange')
ax2.set_title(f"Time Evolution of Electronegativity at Pressure = {p:.2e} torr")
ax2.set_xlabel(r"$\omega_{pi} t$")
ax2.set_ylabel(r"$\alpha$")
ax2.grid(True, which='both', linestyle='--', alpha=0.5)
ax2.legend()
plt.tight_layout()
plt.savefig(pjoin(path_fig, 'time_vs_alpha.pdf'), dpi=1200)

plt.show()
