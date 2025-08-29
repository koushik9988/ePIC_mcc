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
import matplotlib.patches as patches


# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

plot_path = './plots'

path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1



j = NUM_TS
# Plotting
fig, (ax, axden) = plt.subplots(2, 1)

pot = f[f"fielddata/pot/{j}"]
#den_e = f[f"fielddata/den_electron/{j}"]
den_i = f[f"fielddata/den_argon/{j}"]
#rho = den_i[:] - den_e[:]
#den_n = f[f"fielddata/den_negion/{j}"]


x = np.linspace(0, NC, len(pot))
ax.plot(x, pot[:], color='black', label="$\phi$")
#axden.plot(x, den_e[:], color='red', label="$n_{e}$")
axden.plot(x, den_i[:], color='green', label="$n_{e}$")
#axden.plot(x, den_n[:], color='blue', label="$n_{n}$")
ax.legend()
axden.legend()
ax.set_ylabel('$\phi$')
axden.set_ylabel('$density$')
axden.set_xlabel('$x$')
ax.legend(loc='upper right', framealpha=0.5)

# Add ellipse 
sheath_center_x = 508  
sheath_center_y = 0.62 


ellipse_width = 10# Width of the ellipse
ellipse_height = 1.5  # Height of the ellipse

# Add ellipse
ellipse = patches.Ellipse((sheath_center_x, sheath_center_y), width=ellipse_width, height=ellipse_height, 
                          edgecolor='blue', fill=False, linewidth=1, linestyle='--')
axden.add_patch(ellipse)

# Add second ellipse (new one)
sheath_center_x2 = 4.5  # Example new X-position
sheath_center_y2 = 0.62  # Example new Y-position

ellipse_width2 = 10  # Width of the second ellipse
ellipse_height2 = 1.5  # Height of the second ellipse

ellipse2 = patches.Ellipse((sheath_center_x2, sheath_center_y2), width=ellipse_width2, height=ellipse_height2, 
                           edgecolor='black', fill=False, linewidth=1, linestyle='--')
axden.add_patch(ellipse2)


plt.savefig(pjoin(path_fig,'sheath.pdf'),dpi = 1200)

plt.show()
