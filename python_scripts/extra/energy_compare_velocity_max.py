import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py

# hdf5 file name and path
file_name = 'result.h5'

path1 = "../data_run_1"
path2 = "../data_run_2"
path3 = "../data_run_3"
path4 = "../data_run_4"
path5 = "../data_run_5"
path6 = "../data_run_6"

path0 = './plots'

path_fig = pjoin(path1, path0)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# ----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')

# Kinetic energy data
data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]

# Extracting time and potential energy
ts = data1[:, 0]
pe1 = data1[:, 5]
pe2 = data2[:, 5]
pe3 = data3[:, 5]
pe4 = data4[:, 5]
pe5 = data5[:, 5]
pe6 = data6[:, 5]

# Metadata for scaling factors
metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']
metadata_group3 = f3['/metadata']
metadata_group4 = f4['/metadata']
metadata_group5 = f5['/metadata']
metadata_group6 = f6['/metadata']

we1 = metadata_group1.attrs['wpe']
we2 = metadata_group2.attrs['wpe']
we3 = metadata_group3.attrs['wpe']
we4 = metadata_group4.attrs['wpe']
we5 = metadata_group5.attrs['wpe']
we6 = metadata_group6.attrs['wpe']

wp1 = metadata_group1.attrs['wpi']
wp2 = metadata_group2.attrs['wpi']
wp3 = metadata_group3.attrs['wpi']
wp4 = metadata_group4.attrs['wpi']
wp5 = metadata_group5.attrs['wpi']
wp6 = metadata_group6.attrs['wpi']

mfactor1 = wp1/we1
mfactor2 = wp2/we2
mfactor3 = wp3/we3
mfactor4 = wp4/we4
mfactor5 = wp5/we5
mfactor6 = wp6/we6

# Time scaling
ts1 = data1[:, 0] * mfactor1
ts2 = data2[:, 0] * mfactor2
ts3 = data3[:, 0] * mfactor3
ts4 = data4[:, 0] * mfactor4
ts5 = data5[:, 0] * mfactor5
ts6 = data6[:, 0] * mfactor6

# Velocity values
velocities = np.array([5, 10, 15, 20, 25, 30])

# Find the time when PE reaches max for each velocity
times_to_max_pe = [
    ts1[np.argmax(pe1)],
    ts2[np.argmax(pe2)],
    ts3[np.argmax(pe3)],
    ts4[np.argmax(pe4)],
    ts5[np.argmax(pe5)],
    ts6[np.argmax(pe6)]
]

# Plot time to reach max PE as a function of velocity
figsize = np.array([80, 80 / 1.618])  # Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200  # Print resolution
ppi = np.sqrt(1920**2 + 1200**2) / 24  # Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

fig, ax = plt.subplots(figsize=figsize / 10.4, constrained_layout=True, dpi=ppi)

ax.plot(velocities, times_to_max_pe, marker='o', color='blue', label='Time to Max PE')

ax.set_xlabel("$v_{b0}$")
ax.set_ylabel("$\\omega_{{pi}} t$")
ax.grid(True)

ax.legend()
plt.savefig(pjoin(path_fig, 'time_to_max_pe_vs_velocity_max.pdf'), dpi=dpi)
plt.show()
