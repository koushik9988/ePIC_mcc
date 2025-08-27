import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import math as m

# Argument handling
if len(sys.argv) != 5:
    print("Usage: python3 script.py <path1> <path2> <particle_type> <time_step>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
path1 = sys.argv[2]
particle_type = sys.argv[3]
time_step = int(sys.argv[4])

plot_path = './plots/phase'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 files
f1 = h5py.File(pjoin(path, file_name), 'r')
f2 = h5py.File(pjoin(path1, file_name), 'r')
metadata_group1 = f1['/metadata']
metadata_negion1 = f1['/metadata_species/negion']
metadata_negion2 = f2['/metadata_species/negion']


alp1 = metadata_negion1.attrs['density']
alp2 = metadata_negion2.attrs['density']

metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']


eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

wpe1 = metadata_group1.attrs['wpe']
wpi1 = metadata_group1.attrs['wpi']
wpe2 = metadata_group2.attrs['wpe']
wpi2 = metadata_group2.attrs['wpi']

write_interval_phase = metadata_group1.attrs['write_int_phase']
DT_coeff = metadata_group1.attrs['DT_coeff']


mfactor1 = wpi1/wpe1
mfactor2 = wpi2/wpe2


j = time_step * write_interval_phase #(j = TS)

k = j*DT_coeff #wpet
k = k*mfactor1 #wpit

j1 = int(k/(mfactor2*DT_coeff))

j1 = round(j1 / write_interval_phase) * write_interval_phase

print(j)
print(j1)

wpit1 = k 
wpit1 = m.floor(wpit1)
wpit2 = j1*DT_coeff*mfactor2
wpit2 = m.floor(wpit2)



data_phase = f1[f"particle_{particle_type}/{j}"]
datax = data_phase[:, 0]
datavx = data_phase[:, 1]

data_phase1 = f2[f"particle_{particle_type}/{j1}"]
datax1 = data_phase1[:, 0]
datavx1 = data_phase1[:, 1]

# Create the plot
fig, (ax_phase1, ax_phase2) = plt.subplots(2, 1)

# First plot for the first file
#ax_phase1.scatter(datax, datavx, marker='o', color='b', s=1E-1, label=f"$\\alpha = {alp1}$")
ax_phase1.scatter(datax, datavx, marker='o', color='b', s=1, label=f"$\\omega_{{pi}} t = {wpit1:.1f} and \\alpha = {alp1} $")
ax_phase1.set_xlabel('$x$')
ax_phase1.set_ylabel('$v$')
ax_phase1.legend(loc='upper right', framealpha=0.5)


# Second plot for the second file
#ax_phase2.scatter(datax1, datavx1, marker='o', color='b', s=1E-1, label=f"$\\alpha = {alp2}$")
ax_phase2.scatter(datax1, datavx1, marker='o', color='b', s=1, label=f"$\\omega_{{pi}} t = {wpit2:.1f} and \\alpha = {alp2}$")
ax_phase2.set_xlabel('$x$')
ax_phase2.set_ylabel('$v$')
ax_phase2.legend(loc='upper right', framealpha=0.5)
#ax_phase2.set_title(title_text)

# Save the figure
fig_name = f"phase_space_timestep_{wpit1}.png"
plt.savefig(pjoin(path_fig, fig_name),dpi =300)
# Show the plot
plt.show()

