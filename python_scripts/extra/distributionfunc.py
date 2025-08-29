import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os
from os.path import join as pjoin
import sys
import h5py
import seaborn as sns

# Command-line arguments
path = sys.argv[1]
time_step = int(sys.argv[2])

# HDF5 file name and path
file_name = 'result.h5'

plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_negion = f['/metadata_species/negion']

# -----------------------------------------------------------------
eps0 = constants('electric constant')
e = constants('elementary charge')

# Read attributes from hdf5 file 
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
nParticlesE = 500000
nParticlesI = 500000
nParticlesN = 500000
nParticlesB = 500000
Te = 1*e
alp = metadata_negion.attrs['density']
beta = 0.4
n0 = metadata_group.attrs['density']

# Debye length calculation
ne0 = n0 / (1 + alp + beta)
nn0 = alp * ne0
nb0 = beta * ne0

xl = NC
neg_spwt = (nn0 * xl) / nParticlesN
beam_spwt = (nb0 * xl) / nParticlesB

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

# Load data at time step 0
data_phase_n_0 = f["particle_negion/0"]
data_phase_n_0 = data_phase_n_0 
datanvx_0 = data_phase_n_0[:, 1]

data_phase_b_0 = f["particle_beam/0"]
data_phase_b_0 = data_phase_b_0
databvx_0 = data_phase_b_0[:, 1]

# Load data at the specified time step
data_phase_n = f[f"particle_negion/{time_step}"]
data_phase_n = data_phase_n 
datanvx = data_phase_n[:, 1]

data_phase_b = f[f"particle_beam/{time_step}"]
data_phase_b_0 = data_phase_b_0
databvx = data_phase_b[:, 1]

# Plotting the combined distribution function
sns.set_style('whitegrid')

fig, ax_n = plt.subplots(1, 1, figsize=(10, 6))

sns.kdeplot(x=np.concatenate([np.array(datanvx_0), np.array(databvx_0)]), 
            ax=ax_n, color='y', linewidth=2, linestyle='-', fill=True, 
            label=f'Beam + Negative Ions at t=0')
sns.kdeplot(x=np.concatenate([np.array(datanvx), np.array(databvx)]), 
            ax=ax_n, color='y', linewidth=2, linestyle='--', fill=True, 
            label=f'Beam + Negative Ions at t={time_step}')

ax_n.set_xlabel('Velocity (v)')
ax_n.set_ylabel('$f(v)$')
ax_n.legend(loc='lower right', framealpha=0)

plt.tight_layout()
plt.savefig(pjoin(path_fig, f'distributionfunction_t{time_step}.png'))
plt.show()

