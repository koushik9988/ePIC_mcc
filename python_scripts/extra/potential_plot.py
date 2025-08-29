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

# hdf5 file name and path
file_name = 'result.h5'

path = sys.argv[1]


timestep = int(sys.argv[2]) 

plot_path = './plots/field'

path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Read attributes
NC = metadata_group.attrs['NC']

save_fig = metadata_group.attrs['save_fig']

fig,ax = plt.subplots()
j = timestep
pot = f[f"fielddata/pot/{j}"]
ef = f[f"fielddata/efield/{j}"]
x = np.linspace(0, NC, len(pot))
ax.plot(x, pot[:], color='black', label="Potential")
ax.plot(x, ef[:], color='red', label="electric field")
plt.show()