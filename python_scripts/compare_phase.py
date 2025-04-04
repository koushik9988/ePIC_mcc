import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Argument
if len(sys.argv) != 4:
    print("Usage: python3 script.py <path1> <path2> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path1 = sys.argv[1]
path2 = sys.argv[2]
particle_type = sys.argv[3]

plot_path = './plots'

path_fig1 = pjoin(path1, plot_path)
path_fig2 = pjoin(path2, plot_path)

if not os.path.exists(path_fig1):
    os.makedirs(path_fig1)

if not os.path.exists(path_fig2):
    os.makedirs(path_fig2)

# Read hdf5 file
f1 = h5py.File(pjoin(path1, file_name), 'r')
metadata_group1 = f1['/metadata']

f2 = h5py.File(pjoin(path2, file_name), 'r')
metadata_group2 = f2['/metadata']

NC = metadata_group1.attrs['NC']
NUM_TS = metadata_group1.attrs['NUM_TS']
write_interval = metadata_group1.attrs['write_int']
write_interval_phase = metadata_group1.attrs['write_int_phase']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Determine min and max values for phi and v
v_min1, v_max1 = float('inf'), float('-inf')
v_min2, v_max2 = float('inf'), float('-inf')

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    # Get phase space data
    data_phase1 = f1[f"particle_{particle_type}/{j}"]
    datavx1 = data_phase1[:, 1]
    v_min1 = min(v_min1, np.min(datavx1))
    v_max1 = max(v_max1, np.max(datavx1))

    data_phase2 = f2[f"particle_{particle_type}/{j}"]
    datavx2 = data_phase2[:, 1]
    v_min2 = min(v_min2, np.min(datavx2))
    v_max2 = max(v_max2, np.max(datavx2))


# Plotting
fig, (ax_phase1, ax_phase2) = plt.subplots(2, 1)

def animate(i):
    j = i * write_interval_phase

    # Phase space data
    data_phase1 = f1[f"particle_{particle_type}/{j}"]
    datax1 = data_phase1[:, 0]
    datavx1 = data_phase1[:, 1]

    ax_phase1.clear()    
    ax_phase1.scatter(datax1, datavx1, marker='o', color='b', alpha= 1.0, s = 1 , label=f"{particle_type} phase space iaw")
    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.legend(loc='upper right', framealpha=0.5)
    title_text = ' TS: {:.4f}'.format(j)
    ax_phase1.set_title(title_text)

    # Phase space data
    data_phase2 = f2[f"particle_{particle_type}/{j}"]
    datax1 = data_phase2[:, 0]
    datavx1 = data_phase2[:, 1]

    ax_phase2.clear()    
    ax_phase2.scatter(datax1, datavx1, marker='o', color='b', alpha= 1.0, s = 1 , label=f"{particle_type} phase space no iaw")
    ax_phase2.set_xlabel('$x$')
    ax_phase2.set_ylabel('$v$')
    ax_phase2.legend(loc='upper right', framealpha=0.5)
    title_text = ' TS: {:.4f}'.format(j)
    ax_phase1.set_title(title_text)

    return ax_phase1, ax_phase2

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
