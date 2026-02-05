import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# --- Arguments ---
if len(sys.argv) != 4:
    print("Usage: python3 script.py <path> <particle_type> <save_flag>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
save_flag = int(sys.argv[3])  # 1 = save plots, 0 = no save

# --- hdf5 file name and path ---
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# --- Read hdf5 file ---
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# --- Constants and attributes ---
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

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

mFactor = wpi / wpe
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mFactor  # ω_pi * t

# --- Plotting ---
#fig, (ax, ax_phase1) = plt.subplots(2, 1)
fig, ax_phase1 = plt.subplots()

def animate(i):
    j = i * write_interval_phase

    # Phase space data
    data_phase = f[f"particle_{particle_type}/{j}"]
    datax = data_phase[:, 0]
    datavx = data_phase[:, 1]
    ids = data_phase[:, 4].astype(int)  # assuming id is in column 4

    # Color red if collided (id != 0), blue otherwise
    #colors = np.where(ids != 0, 'red', 'blue')

    ax_phase1.clear()

    
    #animate.cbar.remove()
    #animate.cbar = None

    #ax_phase1.scatter(datax, datavx, c=colors, marker='o', s=1, label=rf"{particle_type} phase space ( $\omega_{{pi}}t = {ts[i]:.2f}$)")
    # Normalize colors based on velocity
    vmin, vmax = np.min(datavx), np.max(datavx)
    #sc = ax_phase1.scatter(datax, datavx,c=datavx,cmap='plasma',vmin=vmin,vmax=vmax,s=2, marker='o')
    ax_phase1.scatter(datax, datavx,s=2, marker='o', color='blue')

    # Add colorbar
    #cbar = plt.colorbar(sc, ax=ax_phase1, pad=0.01)
    #cbar.set_label('$v_x$', fontsize=10)
    #cbar.ax.tick_params(labelsize=8)


    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.set_xlim([0, NC])
    #ax_phase1.legend(loc='upper right', framealpha=0.5)

    title_text = f"time : {ts[i]:.2f}, TS: {j}"
    ax_phase1.set_title(title_text)

    # Potential and electric field
    pot = f[f"fielddata/pot/{j}"]
    ef = f[f"fielddata/efield/{j}"]
    x = np.linspace(0, NC, len(pot))

    #ax.clear()
    #ax.plot(x, pot[:], color='red', label=rf"$\phi$ ($\omega_{{pi}}t = {ts[i]:.2f}$)")
    #ax.set_ylabel('$\phi$')
    #ax.legend(loc='upper right', framealpha=0.5)

    # Save if enabled
    if save_flag == 1 and i % 10 == 0:
        fig.savefig(pjoin(path_fig, f"frame_{i:04d}.png"), dpi=150)

    #return ax, ax_phase1
    return ax_phase1

# --- Keyboard controls ---
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
        on_key.ani = animation.FuncAnimation(
            fig, animate, frames=DATA_TS_PHASE,
            blit=False, interval=100, repeat=False
        )
        plt.draw()

on_key.frame = 0
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

animate(on_key.frame)
plt.show()
