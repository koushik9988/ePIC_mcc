import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# --- Arguments ---
if len(sys.argv) != 3:
    print("Usage: python3 animate_space_charge.py <path> <save_flag>")
    sys.exit(1)

path = sys.argv[1]
save_flag = int(sys.argv[2])

# --- Paths ---
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)
os.makedirs(path_fig, exist_ok=True)

# --- Read HDF5 ---
f = h5py.File(pjoin(path, file_name), 'r')
meta = f['/metadata']

# --- Metadata ---
NC = meta.attrs['NC']
NUM_TS = meta.attrs['NUM_TS']
write_interval = meta.attrs['write_int']

DATA_TS = int(NUM_TS / write_interval) + 1
timesteps = [i * write_interval for i in range(DATA_TS)]

# --- Grid ---
num_cells = len(f["fielddata/den_electron/0"])
x = np.linspace(0, NC, num_cells)

# --- Time ---
if "time_var/kinetic_energy" in f:
    ts = f["time_var/kinetic_energy"][:, 0]
else:
    ts = np.arange(DATA_TS)

# ============================================================
#                  FIGURE & AXES (2x2)
# ============================================================
fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

ax_ne, ax_ni = ax[0, 0], ax[0, 1]
ax_rho, ax_both = ax[1, 0], ax[1, 1]

# --- Precompute limits ---
ne_min, ne_max = np.inf, -np.inf
ni_min, ni_max = np.inf, -np.inf
rho_min, rho_max = np.inf, -np.inf

for j in timesteps:
    ne = f[f"fielddata/den_electron/{j}"][:]
    ni = f[f"fielddata/den_ion/{j}"][:]
    rho = ni - ne

    ne_min, ne_max = min(ne_min, ne.min()), max(ne_max, ne.max())
    ni_min, ni_max = min(ni_min, ni.min()), max(ni_max, ni.max())
    rho_min, rho_max = min(rho_min, rho.min()), max(rho_max, rho.max())

# --- Initial plots ---
line_ne, = ax_ne.plot(x, np.zeros_like(x), lw=1.8, color='blue')
line_ni, = ax_ni.plot(x, np.zeros_like(x), lw=1.8, color='red')
line_rho, = ax_rho.plot(x, np.zeros_like(x), lw=1.8, color='purple')

line_ne_both, = ax_both.plot(x, np.zeros_like(x), lw=1.5, color='blue', label=r'$n_e$')
line_ni_both, = ax_both.plot(x, np.zeros_like(x), lw=1.5, color='red',  label=r'$n_i$')

# --- Axis formatting ---
for a in ax.flatten():
    a.set_xlim(0, NC)
    a.grid(True)

ax_ne.set_ylim(1.1 * ne_min, 1.1 * ne_max)
ax_ni.set_ylim(1.1 * ni_min, 1.1 * ni_max)
ax_rho.set_ylim(1.1 * rho_min, 1.1 * rho_max)
ax_both.set_ylim(
    1.1 * min(ne_min, ni_min),
    1.1 * max(ne_max, ni_max)
)

ax_ne.set_ylabel(r'$n_e$')
ax_ni.set_ylabel(r'$n_i$')
ax_rho.set_ylabel(r'$\rho = n_i - n_e$')
ax_rho.set_xlabel(r'$x$')
ax_both.set_xlabel(r'$x$')

ax_ne.set_title("Electron density")
ax_ni.set_title("Ion density")
ax_rho.set_title("Space charge density")
ax_both.set_title(r'$n_e$ and $n_i$')

ax_both.legend(frameon=False)

# ============================================================
#                       ANIMATION
# ============================================================
def animate(i):
    j = timesteps[i]*20

    ne = f[f"fielddata/den_electron/{j}"][:]
    ni = f[f"fielddata/den_ion/{j}"][:]
    rho = ni - ne

    line_ne.set_ydata(ne)
    line_ni.set_ydata(ni)
    line_rho.set_ydata(rho)
    line_ne_both.set_ydata(ne)
    line_ni_both.set_ydata(ni)

    fig.suptitle(f"t = {ts[i]:.3f}   |   TS = {j}", fontsize=12)

    if save_flag == 1 and i % 10 == 0:
        fig.savefig(pjoin(path_fig, f"rho_{i:04d}.png"), dpi=150)

    return (line_ne, line_ni, line_rho, line_ne_both, line_ni_both)

# ============================================================
#                   KEYBOARD CONTROLS
# ============================================================
def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS - 1)
        animate(on_key.frame)
        plt.draw()

    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        animate(on_key.frame)
        plt.draw()

    elif event.key == 'tab':
        on_key.ani = animation.FuncAnimation(
            fig, animate,
            frames=DATA_TS,
            interval=100,
            blit=False,
            repeat=False
        )
        plt.draw()

on_key.frame = 0
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

# --- Initial frame ---
animate(0)
plt.tight_layout()
plt.show()
f.close()
