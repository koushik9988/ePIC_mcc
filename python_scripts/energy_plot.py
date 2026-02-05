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

# Read kinetic energy matrix
data = f["time_var/kinetic_energy"]
ts = data[:, 0]  # First column is time

# Total number of plots: kinetic plots + potential + total
total_plots = spno + 2

# -------------------- Figure 1: KE components + PE + Total Energy --------------------
cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig1, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

total_ke = np.zeros_like(ts)
ke_species_total = {}  # Store total KE of each species for use in Figure 2

# Plot kinetic energy for each species based on species_order
for i, species_name in enumerate(species_order):
    species_group = f["/metadata_species"].get(species_name, None)
    if species_group is None:
        raise ValueError(f"Species '{species_name}' not found in '/metadata_species'.")

    start_idx = 1 + i * 3  # KE_x, KE_y, KE_z
    ke_x = data[:, start_idx]
    ke_y = data[:, start_idx + 1]
    ke_z = data[:, start_idx + 2]

    ke_total = ke_x + ke_y + ke_z
    ke_species_total[species_name] = ke_total

    total_ke += ke_total

    axs[i].plot(ts, ke_x, label=f"${species_name}~KE_x$")
    axs[i].plot(ts, ke_y, label=f"${species_name}~KE_y$")
    axs[i].plot(ts, ke_z, label=f"${species_name}~KE_z$")
    axs[i].set_ylabel(f"{species_name} KE")
    axs[i].legend(loc='upper right', framealpha=0.5)

# Potential energy
pe = data[:, 3 * spno + 1]
potential_energy = pe
axs[spno].plot(ts, pe, label="$pe$")
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno].set_ylabel("Potential Energy")
axs[spno].legend(loc='upper right', framealpha=0.5)

# Total energy
total_energy = total_ke + potential_energy
axs[spno + 1].plot(ts, total_ke, label="KE", color='red')
axs[spno + 1].plot(ts, potential_energy, label="PE", color='green')
axs[spno + 1].plot(ts, total_energy, label="Total Energy", color='black')
axs[spno + 1].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_ylabel("Total Energy")
axs[spno + 1].legend(loc='upper right', framealpha=0.5)

# Hide unused subplots
for i in range(total_plots, len(axs)):
    fig1.delaxes(axs[i])

fig1.suptitle("Figure 1: KE Components, PE and Total Energy", fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig(pjoin(path_fig, "figure1_energy_components.png"))
plt.show()

# -------------------- Figure 2: Total KE per species --------------------
cols2 = int(np.ceil(np.sqrt(spno)))
rows2 = int(np.ceil(spno / cols2))

fig2, axs2 = plt.subplots(rows2, cols2, figsize=(5 * cols2, 3.5 * rows2), sharex=True)
axs2 = axs2.flatten()

for i, species_name in enumerate(species_order):
    ke_total = ke_species_total[species_name]
    axs2[i].plot(ts, ke_total, label=f"${species_name}~KE_{{total}}$", color='blue')
    axs2[i].set_ylabel(f"{species_name} Total KE")
    axs2[i].legend(loc='upper right', framealpha=0.5)

# Hide any unused axes
for i in range(spno, len(axs2)):
    fig2.delaxes(axs2[i])

fig2.suptitle("Figure 2: Total KE per Species", fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig(pjoin(path_fig, "figure2_total_ke_species.png"))
plt.show()

"""
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

# --- Read species order from species metadata --- (as hdf5 can only store alphabetic wise)
species_order = f["/metadata_species"].attrs.get("species_order", None)
if species_order is None:
    raise ValueError("Could not find 'species_order' attribute in /metadata_species.")

print(f"Species order: {species_order}")

# Read kinetic energy matrix
data = f["time_var/kinetic_energy"]
ts = data[:, 0]  # First column is time

# Total number of plots: kinetic plots + potential + total
total_plots = spno + 2

# Calculate grid size (square-like)
cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

total_ke = np.zeros_like(ts)

# Plot kinetic energy for each species based on species_order
for i, species_name in enumerate(species_order):
    # Find index of the species in the file's metadata
    species_group = f["/metadata_species"].get(species_name, None)
    if species_group is None:
        raise ValueError(f"Species '{species_name}' not found in '/metadata_species'.")

    # Extract kinetic energy for the species
    start_idx = 1 + i * 3  # Each species has 3 KE components: x, y, z
    ke_x = data[:, start_idx]
    ke_y = data[:, start_idx + 1]
    ke_z = data[:, start_idx + 2]

    total_ke += ke_x + ke_y + ke_z  # Sum kinetic energies

    axs[i].plot(ts, ke_x, label=f"${species_name}~KE_x$")
    axs[i].plot(ts, ke_y, label=f"${species_name}~KE_y$")
    axs[i].plot(ts, ke_z, label=f"${species_name}~KE_z$")
    axs[i].set_ylabel(f"{species_name} KE")
    axs[i].legend(loc='upper right', framealpha=0.5)

# Potential energy plot
pe = data[:, 3 * spno + 1]
potential_energy = pe

axs[spno].plot(ts, pe, label="$pe$")
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno].set_ylabel("Potential Energy")
axs[spno].legend(loc='upper right', framealpha=0.5)

# Total energy plot
total_energy = total_ke + potential_energy
axs[spno + 1].plot(ts, total_ke, label="KE", color='red')
axs[spno + 1].plot(ts, potential_energy, label="PE", color='green')
axs[spno + 1].plot(ts, total_energy, label="Total Energy", color='black')
#axs[spno + 1].set_ylim(np.min(total_energy) - 5, np.max(total_energy) + 5)
axs[spno + 1].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_ylabel("Total Energy")
axs[spno + 1].legend(loc='upper right', framealpha=0.5)

# Hide unused subplots if any
for i in range(total_plots, len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.show()
"""