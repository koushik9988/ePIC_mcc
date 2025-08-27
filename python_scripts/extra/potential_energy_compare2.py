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
path1 = "../build/data_nocoll"
path2 = "../build/data_1e17"
path3 = "../build/data_5e17"
path4 = "../build/data_1e18"
path5 = "../build/data_1.5e18"
path6 = "../build/data_2e18"
path7 = "../build/data_5e18"

path8 = './plots'

path_fig = pjoin(path6,path5)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')
f7 = h5py.File(pjoin(path7, file_name), 'r')

data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]
data7 = f6["time_var/kinetic_energy"]

ts1 = data1[:,0]
PE1 = data1[:,13]
PE2 = data2[:,13]
PE3 = data3[:,13]
PE4 = data4[:,13]
PE5 = data5[:,13]
PE6 = data6[:,13]
PE7 = data6[:,13]


#---------------plotting -------------------------
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

fig,ax = plt.subplots(figsize= figsize/10.4,constrained_layout=True,dpi=ppi)

ax.plot(ts1, PE1, color ='red',label='no collision')
ax.plot(ts1, PE2, color ='green',label='1e17')
ax.plot(ts1, PE3, color ='blue',label='5e17')
ax.plot(ts1, PE4, color ='purple',label='1e18')
ax.plot(ts1, PE5, color ='violet',label='1.5e18')
ax.plot(ts1, PE6, color = 'black',label='2e18')
#ax.plot(ts1, PE7, color = 'yellow',label='5e18')
ax.set_xlabel('$\omega_{pi}t$')
ax.set_ylabel('Electrostatic potential energy')
ax.set_title("e + H elastic collision")
#ax.set_ylabel('kbe')
ax.grid(True)
ax.legend(loc='upper right',framealpha=0.5)
plt.semilogy()
ax.set_ylim(3e-3,1)
plt.tight_layout()

#if(save_fig == 1):
plt.savefig(pjoin(path1,'pe.pdf'),dpi = dpi)

plt.show()