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
import matplotlib.animation as animation
import scipy.integrate as intg

# hdf5 file name and path 
file_name = 'result.h5'

path1 = "../data_run_1"
path2 = "../data_run_2"
path3 = "../data_run_3"
path4 = "../data_run_4"
path5 = "../data_run_5"
path6 = "../data_run_6"
path7 = "../data_run_7"
path8 = "../data_run_8"
path9 = "../data_run_9"
path0 = './plots'

path_fig = pjoin(path1,path0)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# ----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')
f7 = h5py.File(pjoin(path7, file_name), 'r')
f8 = h5py.File(pjoin(path8, file_name), 'r')
f9 = h5py.File(pjoin(path9, file_name), 'r')

data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]
data7 = f7["time_var/kinetic_energy"]
data8 = f8["time_var/kinetic_energy"]
data9 = f9["time_var/kinetic_energy"]

ts = data1[:,0]
pe1 = data1[:,5]
pe2 = data2[:,5]
pe3 = data3[:,5]
pe4 = data4[:,5]
pe5 = data5[:,5]
pe6 = data6[:,5]
pe7 = data7[:,5]
pe8 = data8[:,5]
pe9 = data9[:,5]


figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

fig, ax = plt.subplots(figsize=figsize/10.4,constrained_layout=True,dpi=ppi)


# Now we create a new plot with alpha on the x-axis and the maximum value of pe on the y-axis
alphas = [0, 0.1, 0.5, 1, 3, 5, 10, 15, 20]
max_pe_values = [np.max(pe1), np.max(pe2), np.max(pe3), np.max(pe4), np.max(pe5), np.max(pe6), np.max(pe7), np.max(pe8), np.max(pe9)]

ax.plot(alphas, max_pe_values, marker='o', linestyle='-', color='b')
ax.set_xlabel("$\\alpha$")#, fontsize=14)
ax.set_ylabel("potential")

#plt.xticks(alphas)  # Set x-ticks to be at each alpha value
plt.grid()

plt.savefig(pjoin(path_fig, 'max_pe_vs_alpha.pdf'), dpi=dpi)
plt.show()
