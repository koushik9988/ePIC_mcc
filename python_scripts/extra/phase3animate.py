import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib as mpp
import matplotlib.animation as animation


# hdf5 file name and paths
file_name = 'result.h5'
path1 = "../data_esw_coll_1e15"
path2 = "../data_esw_coll2_1e15"
path3 = "../data_esw_coll_1e19"

particle_type = sys.argv[1]
#time = int(sys.argv[5])  # Convert time argument to integer

plot_path = './plots'

# Figure setup
figsize = np.array([100, 300 / 1.618])  # Figure size in mm
dpi = 300  # Print resolution
ppi = np.sqrt(1920**2 + 1200**2) / 24  # Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)

# Create subplots
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)

# Read hdf5 files
f1 = h5py.File(pjoin(path1, file_name), 'r')
f2 = h5py.File(pjoin(path2, file_name), 'r')
f3 = h5py.File(pjoin(path3, file_name), 'r')

metadata_group1 = f1['/metadata']
metadata_group2 = f2['/metadata']
metadata_group3 = f3['/metadata']

write_interval_phase = metadata_group1.attrs['write_int_phase']
NUM_TS = metadata_group1.attrs['NUM_TS']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Animation function
def animate(i):
    j = i * write_interval_phase
    ax1.clear()
    ax2.clear()
    ax3.clear()
    
    data_phase1 = f1[f"particle_{particle_type}/{j}"]
    datax1 = data_phase1[:, 0]
    datavx1 = data_phase1[:, 1]
    ax1.scatter(datax1, datavx1, marker='o', color='b', alpha=1.0, s=2, label=f"{particle_type} phase space")
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$v$')
    ax1.legend(loc='upper right', framealpha=0.5)
    ax1.set_title(f'Phase Space from {path1}')

    data_phase2 = f2[f"particle_{particle_type}/{j}"]
    datax2 = data_phase2[:, 0]
    datavx2 = data_phase2[:, 1]
    ax2.scatter(datax2, datavx2, marker='o', color='b', alpha=1.0, s=2, label=f"{particle_type} phase space")
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$v$')
    ax2.legend(loc='upper right', framealpha=0.5)
    ax2.set_title(f'Phase Space from {path2}')

    data_phase3 = f3[f"particle_{particle_type}/{j}"]
    datax3 = data_phase3[:, 0]
    datavx3 = data_phase3[:, 1]
    ax3.scatter(datax3, datavx3, marker='o', color='b', alpha=1.0, s=2, label=f"{particle_type} phase space")
    ax3.set_xlabel('$x$')
    ax3.set_ylabel('$v$')
    ax3.legend(loc='upper right', framealpha=0.5)
    ax3.set_title(f'Phase Space from {path3}')


    return ax1, ax2, ax3

def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS_PHASE - 1)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'tab':
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval= 100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None  # To keep reference to animation

fig.canvas.mpl_connect('key_press_event', on_key)

animate(on_key.frame)
plt.show()
