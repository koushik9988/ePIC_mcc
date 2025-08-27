import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import h5py
from os.path import join as pjoin
import os

# HDF5 file name and path 
file_name = 'result.h5'

# Define paths for data file
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

# Create the output directory if it doesn't exist
path_fig = pjoin(path1, path0)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 files
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')
f4 = h5py.File(pjoin(path4, file_name), 'r')
f5 = h5py.File(pjoin(path5, file_name), 'r')
f6 = h5py.File(pjoin(path6, file_name), 'r')
f7 = h5py.File(pjoin(path7, file_name), 'r')
f8 = h5py.File(pjoin(path8, file_name), 'r')
f9 = h5py.File(pjoin(path9, file_name), 'r')

# Read metadata and calculate scaling factors
we1 = f1['/metadata'].attrs['wpe']
we2 = f2['/metadata'].attrs['wpe']
we3 = f3['/metadata'].attrs['wpe']
we4 = f4['/metadata'].attrs['wpe']
we5 = f5['/metadata'].attrs['wpe']
we6 = f6['/metadata'].attrs['wpe']
we7 = f7['/metadata'].attrs['wpe']
we8 = f8['/metadata'].attrs['wpe']
we9 = f9['/metadata'].attrs['wpe']

wp1 = f1['/metadata'].attrs['wpi']
wp2 = f2['/metadata'].attrs['wpi']
wp3 = f3['/metadata'].attrs['wpi']
wp4 = f4['/metadata'].attrs['wpi']
wp5 = f5['/metadata'].attrs['wpi']
wp6 = f6['/metadata'].attrs['wpi']
wp7 = f7['/metadata'].attrs['wpi']
wp8 = f8['/metadata'].attrs['wpi']
wp9 = f9['/metadata'].attrs['wpi']

mfactor1 = wp1 / we1
mfactor2 = wp2 / we2
mfactor3 = wp3 / we3
mfactor4 = wp4 / we4
mfactor5 = wp5 / we5
mfactor6 = wp6 / we6
mfactor7 = wp7 / we7
mfactor8 = wp8 / we8
mfactor9 = wp9 / we9


data1 = f1["time_var/kinetic_energy"]
data2 = f2["time_var/kinetic_energy"]
data3 = f3["time_var/kinetic_energy"]
data4 = f4["time_var/kinetic_energy"]
data5 = f5["time_var/kinetic_energy"]
data6 = f6["time_var/kinetic_energy"]
data7 = f7["time_var/kinetic_energy"]
data8 = f8["time_var/kinetic_energy"]
data9 = f9["time_var/kinetic_energy"]


pe1 = data1[:, 5]
pe2 = data2[:, 5]
pe3 = data3[:, 5]
pe4 = data4[:, 5]
pe5 = data5[:, 5]
pe6 = data6[:, 5]
pe7 = data7[:, 5]
pe8 = data8[:, 5]
pe9 = data9[:, 5]

#time scaling
max_pe_time1 = data1[np.argmax(pe1), 0] * mfactor1
max_pe_time2 = data2[np.argmax(pe2), 0] * mfactor2
max_pe_time3 = data3[np.argmax(pe3), 0] * mfactor3
max_pe_time4 = data4[np.argmax(pe4), 0] * mfactor4
max_pe_time5 = data5[np.argmax(pe5), 0] * mfactor5
max_pe_time6 = data6[np.argmax(pe6), 0] * mfactor6
max_pe_time7 = data7[np.argmax(pe7), 0] * mfactor7
max_pe_time8 = data8[np.argmax(pe8), 0] * mfactor8
max_pe_time9 = data9[np.argmax(pe9), 0] * mfactor9

# Store the times and corresponding alpha values
max_pe_times = np.array([max_pe_time1, max_pe_time2, max_pe_time3,
                          max_pe_time4, max_pe_time5, max_pe_time6,
                          max_pe_time7, max_pe_time8, max_pe_time9])

# Define alpha values corresponding to each data run
alphas = np.array([0, 0.1, 0.5, 1, 3, 5, 10, 15, 20])


normalized_times = max_pe_times# / np.max(max_pe_times)

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

ax.plot(alphas, normalized_times, marker='o', linestyle='-', color='b')
ax.set_xlabel("$\\alpha$")#, fontsize=14)
ax.set_ylabel("$\\omega_{{pi}} t$")#, fontsize=14)


#plt.xticks(alphas)  # Set x-ticks to be at each alpha value
plt.grid()
plt.savefig(pjoin(path_fig, 'max_pe_time_vs_alpha.pdf'), dpi=dpi)
plt.show()
