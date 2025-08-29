import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Argument
if len(sys.argv) != 3:
    print("Usage: python3 script.py <path1> <path2> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'

path1 = sys.argv[1]
path2 = sys.argv[2]

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
phi_min1, phi_max1 = float('inf'), float('-inf')
phi_min2, phi_max2 = float('inf'), float('-inf')

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    pot1 = f1[f"fielddata/pot/{j}"]
    phi_min1 = min(phi_min1, np.min(pot1))
    phi_max1 = max(phi_max1, np.max(pot1))

    pot2 = f2[f"fielddata/pot/{j}"]
    phi_min2 = min(phi_min2, np.min(pot2))
    phi_max2 = max(phi_max2, np.max(pot2))


# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1)

def animate(i):
    j = i * write_interval_phase

    ax1.clear()
    pot1 = f1[f"fielddata/pot/{j}"]
    ef1 = f1[f"fielddata/efield/{j}"]
    x = np.linspace(0, NC, len(pot1))
    #ax1.plot(x, pot1[:], color='black', label="Potential iaw ")
    ax1.plot(x, ef1[:], color='red', label="efield")
    #ax1.set_ylim([phi_min1,phi_max1])
    ax1.set_ylabel('$\phi$')
    ax1.legend(loc='upper right', framealpha=0.5)
    title_text = ' TS: {:.4f}'.format(j)
    ax1.set_title(title_text)

    ax2.clear()
    pot2 = f2[f"fielddata/pot/{j}"]
    ef2 = f2[f"fielddata/efield/{j}"]
    x = np.linspace(0, NC, len(pot2))
    #ax2.plot(x, pot2[:], color='black', label="Potential no iaw")
    ax2.plot(x, ef[:], color='red', label="efield")
    #ax2.set_ylim([phi_min2,phi_max2])
    ax2.set_ylabel('$\phi$')
    ax2.legend(loc='upper right', framealpha=0.5)
    
    return ax1 , ax2

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
