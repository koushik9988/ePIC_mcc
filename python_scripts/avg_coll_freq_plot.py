import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

file_name = "result.h5"
path = sys.argv[1]
output_dir = "./plots"
path_fig = pjoin(path, output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), "r")

# --- Read number of species from metadata ---
spno = f["/metadata"].attrs.get("spno", None)
species_order = f["/metadata_species"].attrs.get("species_order", None)
print(f"Species order: {species_order}")
#Read average collision frequency dataset ---
data = f["time_var/avg_collision_frequency"] 
ts = data[:, 0]  # first column = time

# Plot average collision frequency per species
cols = int(np.ceil(np.sqrt(spno-1)))
rows = int(np.ceil((spno-1) / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

#for i, species_name in enumerate(species_order):
for i, species_name in enumerate(species_order[:-1]):
    avg_freq = data[:, i + 1]  # column i+1 corresponds to species i
    axs[i].plot(ts, avg_freq, label=f"{species_name} avg coll freq", color="blue")
    axs[i].set_ylabel("⟨$n_g \sigma v$⟩")
    axs[i].legend(loc="upper right", framealpha=0.5)

# Hide unused axes
for i in range(spno-1, len(axs)):
    fig.delaxes(axs[i])

fig.suptitle("Average Collision Frequency per Species", fontsize=16)
plt.tight_layout()
plt.savefig(pjoin(path_fig, "avg_collision_freq.png"))
plt.show()
