import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib as mpp

# hdf5 file name and paths
file_name = 'result.h5'

path1 = "../data_esw_coll_1e17"
path2 = "../data_esw_coll_1e19"
path3 = "../data_esw_nocoll"

particle_type = sys.argv[1]
time = int(sys.argv[2])  # Convert time argument to integer

plot_path = './plots'
path_fig1 = pjoin(path1, plot_path)
path_fig2 = pjoin(path2, plot_path)
path_fig3 = pjoin(path3, plot_path)

if not os.path.exists(path_fig1):
    os.makedirs(path_fig1)
if not os.path.exists(path_fig2):
    os.makedirs(path_fig2)
if not os.path.exists(path_fig3):
    os.makedirs(path_fig3)

# Read hdf5 files
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')

metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']
metadata_group3 = f3['/metadata']


LD = metadata_group1.attrs['LDe']
LDi = metadata_group1.attrs['LDi']
wpe = metadata_group1.attrs['wpe']
wpi = metadata_group1.attrs['wpi']
DT_coeff = metadata_group1.attrs['DT_coeff']

write_interval_phase = metadata_group1.attrs['write_int_phase']

# Compute time index
j = time * write_interval_phase

wpeit = j * DT_coeff

#mFactor = wpi/wpe
#wpit = mFactor * wpeit

# Figure setup
#figsize = np.array([100, 300 / 1.618])  # Figure size in mm
dpi = 1200  # Print resolution
ppi = np.sqrt(1920**2 + 1200**2) / 24  # Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)

# Create subplots
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1,figsize=(15, 8))

# Load phase space data and plot for path1
data_phase1 = f1[f"particle_{particle_type}/{j}"]
datax1 = data_phase1[:, 0]
datavx1 = data_phase1[:, 1]
ax1.scatter(datax1, datavx1, marker='o', color='b', alpha=1.0, s=1, label=f"{particle_type} phase space (1e17)")
ax1.set_xlabel('$x$')
ax1.set_ylabel('$v$')
title_text = 'wpet : {:.2f}\n TS: {:.4f}'.format(wpeit,j)
ax1.set_title(title_text)
ax1.legend(loc='upper right', framealpha=0.5)


# Load phase space data and plot for path2
data_phase2 = f2[f"particle_{particle_type}/{j}"]
datax2 = data_phase2[:, 0]
datavx2 = data_phase2[:, 1]
ax2.scatter(datax2, datavx2, marker='o', color='b', alpha=1.0, s=1, label=f"{particle_type} phase space (1e19)")
ax2.set_xlabel('$x$')
ax2.set_ylabel('$v$')
ax2.legend(loc='upper right', framealpha=0.5)


# Load phase space data and plot for path3
data_phase3 = f3[f"particle_{particle_type}/{j}"]
datax3 = data_phase3[:, 0]
datavx3 = data_phase3[:, 1]
ax3.scatter(datax3, datavx3, marker='o', color='b', alpha=1.0, s=1, label=f"{particle_type} phase space (no collision)")
ax3.set_xlabel('$x$')
ax3.set_ylabel('$v$')
ax3.legend(loc='upper right', framealpha=0.5)


# Save and show the plot
plt.savefig(pjoin(path1, plot_path, f'phasespace_{particle_type}_t{time}.png'), dpi=dpi)
plt.show()

# Close files
f1.close()
f2.close()
f3.close()
