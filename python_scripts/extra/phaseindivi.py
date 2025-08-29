import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation
import matplotlib as mpp

# Argument
if len(sys.argv) != 4:
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
particle_type1 = sys.argv[2]
particle_type2 = sys.argv[3]
#time = sys.argv[3]

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

#j = time*write_interval_phase

#wpeit = j*DT_coeff 
#wpit = mFactor*wpeit

figsize = np.array([100,100/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)



# Plotting
fig , ax_phase1 = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
"""
data_phase = f[f"particle_{particle_type1}/{j}"]
datax = data_phase[:, 0]
datavx = data_phase[:, 1]


data_phase = f[f"particle_{particle_type2}/{j}"]
datax = data_phase[:, 0]
datavx = data_phase[:, 1]


den = f[f"fielddata/den_{particle_type1}/{j}"]
x = np.linspace(0, NC, len(den))


ax_phase1.scatter(datax, datavx, marker='.', color='b', alpha=1.0, s = 1 , label=f"{particle_type} phase space")
ax_phase1.set_xlabel('$x$')
ax_phase1.set_ylabel('$v$')
ax_phase1.legend(loc='upper right', framealpha=0.5)

plt.savefig(pjoin(path,'phasespace.pdf'),dpi=dpi)   
plt.show()

"""
def animate(i):
    j = i * write_interval_phase

    # Phase space data
    data_phase = f[f"particle_{particle_type1}/{j}"]
    datax = data_phase[:, 0]
    datavx = data_phase[:, 1]

    data_phase1 = f[f"particle_{particle_type2}/{j}"]
    datax1 = data_phase1[:, 0]
    datavx1 = data_phase1[:, 1]

    #den = f[f"fielddata/den_{particle_type1}/{j}"]
    #x = np.linspace(0, NC, len(den))

    ax_phase1.clear()
    
    ax_phase1.scatter(datax, datavx, marker='o', color='b', alpha= 1.0, s = 2 , label=f"{particle_type1} phase space")
    ax_phase1.scatter(datax1, datavx1, marker='o', color='r', alpha= 1.0, s = 2 , label=f" {particle_type2} phase space")
    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.legend(loc='upper right', framealpha=0.5)


    return ax_phase1

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
