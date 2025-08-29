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


#hdf5 file name and path 
file_name = 'result.h5'

path1 = "../data_run_1"
path7 = "../data_run_7"
path8 = "../data_run_8"
path9 = "../data_run_9"
path0 = './plots'


path_fig = pjoin(path1,path0)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f7 = h5py.File(pjoin(path7, file_name), 'r')
f8 = h5py.File(pjoin(path8, file_name), 'r')
f9 = h5py.File(pjoin(path9, file_name), 'r')

metadata_group1 = f1['/metadata']
metadata_group7 = f7['/metadata']
metadata_group8 = f8['/metadata']
metadata_group9 = f9['/metadata']

we1 = metadata_group1.attrs['wpe']
we7 = metadata_group7.attrs['wpe']
we8 = metadata_group8.attrs['wpe']
we9 = metadata_group9.attrs['wpe']

wp1 = metadata_group1.attrs['wpi']
wp7 = metadata_group7.attrs['wpi']
wp8 = metadata_group8.attrs['wpi']
wp9 = metadata_group9.attrs['wpi']

mfactor1 = wp1/we1
mfactor7 = wp7/we7
mfactor8 = wp8/we8
mfactor9 = wp9/we9

data1 = f1["time_var/kinetic_energy"]
data7 = f7["time_var/kinetic_energy"]
data8 = f8["time_var/kinetic_energy"]
data9 = f9["time_var/kinetic_energy"]

ts1 = data1[:,0]*mfactor1
ts7 = data7[:,0]*mfactor7
ts8 = data8[:,0]*mfactor8
ts9 = data9[:,0]*mfactor9


pe1 = data1[:,5]
pe7 = data7[:,5]
pe8 = data8[:,5]
pe9 = data9[:,5]

#ts *= mFactor
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

ax.plot(ts1, pe1, label="$\\beta = 0.4$", color='black') 
ax.plot(ts7, pe7, label="$\\beta = 0.6$", color='green')
ax.plot(ts8, pe8, label="$\\beta = 0.8$", color='red')
ax.plot(ts9, pe9, label="$\\beta = 1.0$", color='blue')
ax.set_xlabel("$\\omega_{{pi}} t$",fontsize=14)
ax.set_ylabel(r"$\frac{\frac{1}{2} \epsilon_0 E^2}{Th_{ion}}$",fontsize=14)

plt.grid()
ax.legend()
plt.savefig(pjoin(path_fig,'pe_comp_beta.pdf'),dpi = dpi)
plt.show()
