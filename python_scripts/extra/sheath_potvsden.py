import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

start_idx = int(sys.argv[2])



plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
DT_coeff = metadata_group.attrs['DT_coeff']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
j = NUM_TS

fig, ax = plt.subplots(figsize=(8, 6))

pot = f[f"fielddata/pot/{j}"]
den_e = f[f"fielddata/den_electron/{j}"]
den_i = f[f"fielddata/den_ion/{j}"]
#den_n = f[f"fielddata/den_negion/{j}"]

end_idx = NC


# slice the arrays
x = np.linspace(0, NC, len(pot))  
pot_range = pot[start_idx:end_idx+1] 
den_e_range = den_e[start_idx:end_idx+1] 
den_i_range = den_i[start_idx:end_idx+1]  
#den_n_range = den_n[start_idx:end_idx+1]  

# Make potential at start_idx zero
pot_ref = pot_range[0]  
pot_sub = pot_range - pot_ref  


ax.plot(pot_sub, den_e_range, label="$n_{e}$", color='red', linestyle='-', linewidth=2)
ax.plot(pot_sub, den_i_range, label="$n_{i}$", color='green', linestyle='-', linewidth=2)
#ax.plot(pot_sub, den_n_range, label="$n_{n}$", color='black', linestyle='-', linewidth=2)

ax.set_xlabel(f'$\phi$')
ax.set_ylabel('$density$', fontsize=10)
ax.set_xlim(ax.get_xlim()[::-1])  # Reverse the x-axis
ax.legend(loc='upper right', framealpha=0.5)

plt.savefig(pjoin(path_fig,'sheathdenvspot.pdf'),dpi = 1200)

#plt.tight_layout()
plt.show()
