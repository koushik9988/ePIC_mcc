import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys


if len(sys.argv) != 4:
    print("Usage: python3 plot_xt_or_vxt.py <path> <particle_type> <plot_type>")
    print("plot_type: x-t, vx-t, vy-t, vz-t")
    sys.exit(1)

#  Read command-line arguments
path = sys.argv[1]
particle_type = sys.argv[2]
plot_type = sys.argv[3]

file_name = 'result.h5'
file_path = os.path.join(path, file_name)
plot_dir = os.path.join(path, 'plots')
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

#Open HDF5 file
f = h5py.File(file_path, 'r')

#Read metadata
metadata = f['/metadata'].attrs
NC = metadata['NC']
NUM_TS = metadata['NUM_TS']
write_interval = metadata['write_int_phase']
DT = metadata['DT_coeff']


steps = int(NUM_TS / write_interval) + 1

# Set column index for selected plot_type 
if plot_type == 'x-t':
    col_index = 0
elif plot_type == 'vx-t':
    col_index = 1
elif plot_type == 'vy-t':
    col_index = 2
elif plot_type == 'vz-t':
    col_index = 3
else:
    print("Invalid plot_type! Use one of: x-t, vx-t, vy-t, vz-t")
    sys.exit(1)

#  Gather all data for value range 
all_vals = []
for i in range(steps):
    time_step = i * write_interval
    data = f[f"particle_{particle_type}/{time_step}"]
    vals = data[:, col_index]
    all_vals.append(vals)

# Combine all values 
all_vals = np.concatenate(all_vals)


#vmin = np.percentile(all_vals, 1)
#vmax = np.percentile(all_vals, 99)

vmin = np.min(all_vals)
vmax = np.max(all_vals)



#  Create bins for y-axis (x or velocity) 
nbins_y = 200
ybins = np.linspace(vmin, vmax, nbins_y + 1)

#  Create 2D histogram: y-axis vs time 
hist2d = np.zeros((nbins_y, steps), dtype=np.int32)

for i in range(steps):
    time_step = i * write_interval
    data = f[f"particle_{particle_type}/{time_step}"]
    vals = data[:, col_index]
    hist, _ = np.histogram(vals, bins=ybins)
    hist2d[:, i] = hist


f.close()


fig, ax = plt.subplots(figsize=(8, 6))
time_max = NUM_TS * DT
extent = [0, time_max, ybins[0], ybins[-1]]

if plot_type == 'x-t':
    y_label = '$x$'
elif plot_type == 'vx-t':
    y_label = '$v_x$'
elif plot_type == 'vy-t':
    y_label = '$v_y$'
elif plot_type == 'vz-t':
    y_label = '$v_z$'

img = ax.imshow(hist2d, extent=extent, origin='lower', interpolation= 'bilinear', aspect='auto', cmap='hot')
ax.set_xlabel('Time [$\omega_{pe}^{-1}$]')
ax.set_ylabel(y_label)
#ax.set_title(f"{y_label} vs Time for {particle_type}")

cbar = plt.colorbar(img, ax=ax)
if plot_type == 'x-t':
    cbar.set_label("$\int f(x,v) dx$")
else:
    cbar.set_label("$\int f(x,v) dv_x$")



output_name = f"{plot_type}_{particle_type}.png"
output_path = os.path.join(plot_dir, output_name)
plt.savefig(output_path, dpi=300)
print(f"Saved plot to {output_path}")


plt.show()
