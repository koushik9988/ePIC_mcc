import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# --- Command-line arguments ---
if len(sys.argv) < 4:
    print("Usage: python script.py <path> <species1> <species2> [--save]")
    sys.exit(1)

path = sys.argv[1]
species1 = sys.argv[2]
species2 = sys.argv[3]
save_flag = int(sys.argv[4])  # Check for save flag

# --- Paths and HDF5 file ---
file_name = 'result.h5'
plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# --- Constants ---
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# --- Simulation parameters ---
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
normscheme = metadata_group.attrs['norm_scheme']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:", mfactor)

# --- Time axis ---
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor  # omega_pi * t

# --- Get species weights ---
metadata_s1 = f.get(f'/metadata_species/{species1}')
metadata_s2 = f.get(f'/metadata_species/{species2}')
if metadata_s1 is None or metadata_s2 is None:
    print(f"Metadata for '{species1}' or '{species2}' not found.")
    sys.exit(1)

spwt1 = metadata_s1.attrs['spwt']
spwt2 = metadata_s2.attrs['spwt']

# --- Determine velocity and spatial bounds ---
v_min, v_max = float('inf'), float('-inf')
x_min, x_max = float('inf'), float('-inf')
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    v1 = f[f"particle_{species1}/{j}"][:, 1]
    v2 = f[f"particle_{species2}/{j}"][:, 1]
    x1 = f[f"particle_{species1}/{j}"][:, 0]
    x2 = f[f"particle_{species2}/{j}"][:, 0]
    combined_v = np.concatenate([v1, v2])
    combined_x = np.concatenate([x1, x2])
    v_min = min(v_min, np.min(combined_v))
    v_max = max(v_max, np.max(combined_v))
    x_min = min(x_min, np.min(combined_x))
    x_max = max(x_max, np.max(combined_x))

# --- Plot setup ---
#fig, (ax_ps, ax_v) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'hspace': 0.1})
fig, (ax_ps, ax_v) = plt.subplots(2, 1, figsize=(8, 10), sharex=True, gridspec_kw={'hspace': 0.1, 'height_ratios': [3, 1]})

# --- Initial data (t=0) ---
j0 = 0
x1_0 = f[f"particle_{species1}/{j0}"][:, 0]
v1_0 = f[f"particle_{species1}/{j0}"][:, 1]
x2_0 = f[f"particle_{species2}/{j0}"][:, 0]
v2_0 = f[f"particle_{species2}/{j0}"][:, 1]

# For f(v)
nbins_v_high = 500
hist1_v0, bins_v = np.histogram(v1_0, bins=nbins_v_high, range=(v_min, v_max), density=True)
hist2_v0, _ = np.histogram(v2_0, bins=nbins_v_high, range=(v_min, v_max), density=True)
hist1_v0 *= spwt1
hist2_v0 *= spwt2
hist_combined_v0 = (hist1_v0 + hist2_v0) / 2
hist_combined_v0 = hist_combined_v0 / np.max(hist_combined_v0) if np.max(hist_combined_v0) > 0 else hist_combined_v0
bin_centers_v0 = (bins_v[:-1] + bins_v[1:]) / 2

# --- Animation function ---
def animate(i):
    j = i * write_interval_phase
    x1 = f[f"particle_{species1}/{j}"][:, 0]
    v1 = f[f"particle_{species1}/{j}"][:, 1]
    x2 = f[f"particle_{species2}/{j}"][:, 0]
    v2 = f[f"particle_{species2}/{j}"][:, 1]

    # Velocity distribution f(v)
    nbins_v_low = 100
    hist1_v, bins_v = np.histogram(v1, bins=nbins_v_low, range=(v_min, v_max), density=True)
    hist2_v, _ = np.histogram(v2, bins=nbins_v_low, range=(v_min, v_max), density=True)
    hist1_v *= spwt1
    hist2_v *= spwt2
    hist_combined_v = (hist1_v + hist2_v) / 2
    hist_combined_v = hist_combined_v / np.max(hist_combined_v) if np.max(hist_combined_v) > 0 else hist_combined_v
    bin_centers_v = (bins_v[:-1] + bins_v[1:]) / 2

    # Update phase space scatter plot (v on x-axis, x on y-axis)
    ax_ps.clear()
    ax_ps.scatter(v1, x1, s=1, color='blue', alpha=0.3, label="electron")
    ax_ps.scatter(v2, x2, s=1, color='red', alpha=0.3, label="electron beam")
    ax_ps.set_ylabel('$x$')
    #ax_ps.set_xlabel('$v$')
    ax_ps.set_ylim(x_min, x_max)
    #ax_ps.set_xlim(v_min, v_max)
    ax_ps.set_xlim(-6, 11)
    #ax_ps.set_title(f'Phase Space at $\omega_{{pe}}t = {ts[i]:.2f}$')
    ax_ps.legend(loc='upper right', framealpha=0.5)
    ax_ps.grid(True, alpha=0.3)

    # Update f(v) plot
    ax_v.clear()
    ax_v.plot(bin_centers_v0, hist_combined_v0, linestyle='--', color='gray', label=r"$\omega_{pe}t = 0$")
    ax_v.plot(bin_centers_v, hist_combined_v, color='blue', label=rf"$\omega_{{pe}}t = {ts[i]:.2f}$")
    ax_v.set_xlabel('$v$')
    ax_v.set_ylabel('$f(v)$')
    #ax_v.set_xlim(v_min, v_max)
    ax_v.set_xlim(-6, 11)
    ax_v.set_ylim(0, 1.1)
    #ax_v.set_title('Velocity Distribution')
    ax_v.legend(loc='upper right', framealpha=0.5)
    ax_v.grid(True, alpha=0.3)

    # Save if flag enabled
    if save_flag==1 and i % 10 == 0:
        fig.savefig(pjoin(path_fig, f"frame_{i:04d}.png"), dpi=500, bbox_inches='tight')

    return ax_ps, ax_v

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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

# --- Start ---
animate(on_key.frame)
plt.show()