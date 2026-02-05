import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

# --- Setup ---
if len(sys.argv) != 3:
    raise ValueError("Usage: python compare_energy.py <path1> <path2>")

paths = [sys.argv[1], sys.argv[2]]
labels = ["Run 1", "Run 2"]
styles = ['-', '--']  # Line styles for runs

file_name = 'result.h5'
output_dir = './plots'
path_fig = pjoin(paths[0], output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# --- Open both files ---
files = [h5py.File(pjoin(p, file_name), 'r') for p in paths]

# --- Metadata ---
spno = files[0]["/metadata"].attrs.get("spno", None)
species_order = files[0]["/metadata_species"].attrs.get("species_order", None)

if spno is None or species_order is None:
    raise ValueError("Missing 'spno' or 'species_order' attribute in metadata.")

print(f"Number of species : {spno}")
print(f"Species order     : {species_order}")

# --- Setup plotting grid ---
total_plots = spno + 2
cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

# Define colors and line styles for each run
styles = ['-', '--']
colors = ['tab:blue', 'tab:orange']  # One color per run

for run_idx, f in enumerate(files):
    data = f["time_var/kinetic_energy"]
    ts = data[:, 0]  # Time array

    total_ke = np.zeros_like(ts)

    # --- Kinetic Energy Plots ---
    for sp_idx, species_name in enumerate(species_order):
        ke_x = data[:, 1 + sp_idx * 3]
        ke_y = data[:, 2 + sp_idx * 3]
        ke_z = data[:, 3 + sp_idx * 3]

        total_ke += ke_x + ke_y + ke_z

        axs[sp_idx].plot(ts, ke_x, linestyle=styles[run_idx], color=colors[run_idx], label=f"{labels[run_idx]}: KE_x")
        axs[sp_idx].plot(ts, ke_y, linestyle='dotted', color=colors[run_idx], label=f"{labels[run_idx]}: KE_y")
        axs[sp_idx].plot(ts, ke_z, linestyle='dashdot', color=colors[run_idx], label=f"{labels[run_idx]}: KE_z")

        axs[sp_idx].set_ylabel(f"{species_name} KE")
        axs[sp_idx].legend(loc='upper right', framealpha=0.5)

    # --- Potential Energy ---
    pex = data[:, 3 * spno + 1]
    axs[spno].plot(ts, pex, styles[run_idx], color=colors[run_idx], label=f"{labels[run_idx]}: $\int E_x^2$")
    axs[spno].set_ylabel("Potential Energy")
    axs[spno].legend(loc='upper right', framealpha=0.5)

    # --- Total Energy ---
    total_energy = total_ke + pex
    axs[spno + 1].plot(ts, total_energy, styles[run_idx], color=colors[run_idx], label=f"{labels[run_idx]}: Total")
    axs[spno + 1].set_ylim(np.min(total_energy) - 5, np.max(total_energy) + 5)
    axs[spno + 1].set_ylabel("Total Energy")
    axs[spno + 1].legend(loc='upper right', framealpha=0.5)


# --- Label x-axis for bottom plots ---
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_xlabel("$\omega_{pe}t$")

# --- Clean extra subplots ---
for i in range(total_plots, len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.show()
