import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys

# Usage check
if len(sys.argv) != 4:
    print("Usage: python3 plot_xt_or_vxt.py <path> <particle_type> <plot_type>")
    print("plot_type options: x-t, vx-t, vy-t, vz-t")
    sys.exit(1)

# Args
file_name = 'result.h5'
path = sys.argv[1]
particle_type = sys.argv[2]
plot_type = sys.argv[3]

plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read file and metadata
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

nbins_y = 200 
nbins_t = DATA_TS_PHASE

# Select data column
column_map = {
    'x-t': 0,
    'vx-t': 1,
    'vy-t': 2,
    'vz-t': 3,
}
if plot_type not in column_map:
    print("Invalid plot_type! Use one of: 'x-t', 'vx-t', 'vy-t', 'vz-t'")
    sys.exit(1)

col_idx = column_map[plot_type]

# Collect all values to determine histogram range
all_vals = []
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data = f[f"particle_{particle_type}/{j}"]
    all_vals.append(data[:, col_idx])
all_vals = np.concatenate(all_vals)

# Set histogram bounds
vmin, vmax = np.percentile(all_vals, [1, 99])
y_bins = np.linspace(vmin, vmax, nbins_y + 1)

# Create histogram
hist2d = np.zeros((nbins_y, nbins_t), dtype=np.int32)
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data = f[f"particle_{particle_type}/{j}"]
    vals = data[:, col_idx]
    hist, _ = np.histogram(vals, bins=y_bins)
    hist2d[:, i] = hist

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))
extent = [0, NUM_TS * DT_coeff, y_bins[0], y_bins[-1]]
aspect = 'auto'
ylabel_map = {
    'x-t': '$x$',
    'vx-t': '$v_x$',
    'vy-t': '$v_y$',
    'vz-t': '$v_z$',
}
ylabel = ylabel_map[plot_type]

im = ax.imshow(hist2d, extent=extent, origin='upper', interpolation='bilinear', aspect=aspect, cmap='hot')
ax.set_xlabel(r'Time [$\omega_{pe}^{-1}$]')
ax.set_ylabel(ylabel)
ax.set_title(f"{ylabel} vs Time for {particle_type}")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Particle Count per Bin")

# Save and show
outname = f"{plot_type}_{particle_type}.png"
plt.savefig(pjoin(path_fig, outname), dpi=300)
print(f"Saved plot to {pjoin(path_fig, outname)}")
plt.show()
