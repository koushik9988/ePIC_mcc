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
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
particle_type = sys.argv[2]

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
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
wpe = metadata_group.attrs['wpe']
wpi = metadata_group.attrs['wpi']
save_fig = metadata_group.attrs['save_fig']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mFactor = wpi/wpe
data = f["time_var/kinetic_energy"]
ts = data[:,0]
ts *= mFactor # converted to w_pi*t

# Determine min and max values for phi and v
phi_min, phi_max = float('inf'), float('-inf')
v_min, v_max = float('inf'), float('-inf')

"""
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    # Get potential data
    pot = f[f"fielddata/pot/{j}"]
    phi_min = min(phi_min, np.min(pot))
    phi_max = max(phi_max, np.max(pot))
    pot = f[f"fielddata/efield/{j}"]
    ef_min = min(phi_min, np.min(pot))
    ef_max = max(phi_max, np.max(pot))

    # Get phase space data
    data_phase = f[f"particle_{particle_type}/{j}"]
    datavx = data_phase[:, 1]
    v_min = min(v_min, np.min(datavx))
    v_max = max(v_max, np.max(datavx))
"""

# Plotting
fig, (ax, ax_phase1) = plt.subplots(2, 1)

def animate(i):
    j = i * write_interval_phase

    # Phase space data
    data_phase = f[f"particle_{particle_type}/{j}"]
    datax = data_phase[:, 0]
    datavx = data_phase[:, 1]
    ids = data_phase[:, 4].astype(int)  # assuming id is in column 4

    # Example: color red if collided (id != 0), blue if not
    colors = np.where(ids != 0, 'red', 'blue')

    den = f[f"fielddata/den_{particle_type}/{j}"]
    x = np.linspace(0, NC, len(den))

    ax_phase1.clear()
    ax_phase1.scatter(datax, datavx, c=colors, marker='o', s=1,
                      label=f"{particle_type} phase space")
    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.set_xlim([0, NC])
    ax_phase1.legend(loc='upper right', framealpha=0.5)
    title_text = 'time : {:.2f}\n TS: {:.4f}'.format(ts[i], j)
    ax_phase1.set_title(title_text)

    ax.clear()
    pot = f[f"fielddata/pot/{j}"]
    ef = f[f"fielddata/efield/{j}"]
    x = np.linspace(0, NC, len(pot))
    ax.plot(x, pot[:], color='red', label="efield")
    ax.set_ylabel('$\phi$')
    ax.legend(loc='upper right', framealpha=0.5)

    return ax, ax_phase1


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
