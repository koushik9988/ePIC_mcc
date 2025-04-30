import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

file_name = 'result.h5'
path = sys.argv[1]
output_dir = './plots'
path_fig = pjoin(path, output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')

# --- Read number of species from species metadata section ---
spno = f["/metadata"].attrs.get("spno", None)
if spno is None:
    raise ValueError("Could not find 'spno' attribute in /metadata.")

print(f"Number of species : {spno}")

# --- Read species order from species metadata ---
species_order = f["/metadata_species"].attrs.get("species_order", None)
if species_order is None:
    raise ValueError("Could not find 'species_order' attribute in /metadata_species.")

print(f"Species order: {species_order}")

# Read kinetic energy 
data = f["time_var/kinetic_energy"]
ts = data[:, 0]  # First column is time

# Total number of plots: kinetic plots + potential + total
total_plots = spno + 2


cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

total_ke = np.zeros_like(ts)


for i, species_name in enumerate(species_order):
    ke_idx = 1 + i  
    ke = data[:, ke_idx]
    total_ke += ke

    axs[i].plot(ts, ke, label=f"{species_name} KE")
    axs[i].set_ylabel(f"{species_name} KE")
    axs[i].legend(loc='upper right', framealpha=0.5)

pe_idx = 1 + spno 
pe = data[:, pe_idx]

axs[spno].plot(ts, pe, label="Potential Energy")
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno].set_ylabel("Potential Energy")
axs[spno].legend(loc='upper right', framealpha=0.5)

# Total energy plot
total_energy = total_ke + pe
axs[spno + 1].plot(ts, total_energy, label="Total Energy", color='black')
axs[spno + 1].set_ylim(np.min(total_energy) - 5, np.max(total_energy) + 5)
axs[spno + 1].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_ylabel("Total Energy")
axs[spno + 1].legend(loc='upper right', framealpha=0.5)

# Hide unused subplots
for i in range(total_plots, len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.savefig(pjoin(path, 'ke_pe.pdf'), dpi=900)
plt.show()
