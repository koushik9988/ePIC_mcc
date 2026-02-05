import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin

# --------------------- Argument Handling ---------------------
if len(sys.argv) != 2:
    print("Usage: python3 entropy_vs_time_all_species.py <data_path>")
    sys.exit(1)

path = sys.argv[1]
file_name = 'result.h5'
full_path = pjoin(path, file_name)

# --------------------- Create Output Directory ---------------------
output_dir = './plots'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

# --------------------- Open HDF5 File ---------------------
f = h5py.File(full_path, 'r')

# --------------------- Simulation Metadata ---------------------
metadata = f["/metadata"]
spno = metadata.attrs.get("spno", None)
if spno is None:
    raise ValueError("Missing 'spno' attribute in /metadata.")

write_interval_phase = metadata.attrs["write_int_phase"]
NUM_TS = metadata.attrs["NUM_TS"]
wpe = metadata.attrs["wpe"]
wpi = metadata.attrs["wpi"]
normscheme = metadata.attrs["norm_scheme"]
mfactor = wpi / wpe if normscheme == 0 else 1.0
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1
ts = f["time_var/kinetic_energy"][:, 0] * mfactor

# --------------------- Species Order ---------------------
species_order = f["/metadata_species"].attrs.get("species_order", None)
if species_order is None:
    raise ValueError("Missing 'species_order' in /metadata_species.")

# --------------------- Entropy Function ---------------------
def compute_entropy(data, nbins, vmin, vmax, spwt):
    hist, _ = np.histogram(data, bins=nbins, range=(vmin, vmax), density=True)
    hist *= spwt
    hist = np.where(hist > 0, hist, 1e-20)  # avoid log(0)
    dv = (vmax - vmin) / nbins
    return -np.sum(hist * np.log(hist)) * dv

# --------------------- Plotting Setup ---------------------
plt.figure(figsize=(10, 6))
nbins = 100

# --------------------- Total Entropy Initialization ---------------------
S_total_t0 = 0.0
S_total_t = [0.0 for _ in range(DATA_TS_PHASE)]  # will hold total entropy at each t

# --------------------- Loop Over All Species ---------------------
for species_name in species_order:
    print(f"Processing species: {species_name}")
    species_group = f[f"/metadata_species/{species_name}"]
    spwt = species_group.attrs["spwt"]

    # Determine v_min and v_max across all frames for this species
    v_min, v_max = float('inf'), float('-inf')
    for i in range(DATA_TS_PHASE):
        j = i * write_interval_phase
        v = f[f"particle_{species_name}/{j}"][:, 1]
        v_min = min(v_min, np.min(v))
        v_max = max(v_max, np.max(v))

    # Compute S0
    data0 = f[f"particle_{species_name}/0"][:, 1]
    S0 = compute_entropy(data0, nbins, v_min, v_max, spwt)
    S_total_t0 += S0  # accumulate S(0)

    # Compute S(t)
    delta_S = []
    for i in range(DATA_TS_PHASE):
        j = i * write_interval_phase
        v = f[f"particle_{species_name}/{j}"][:, 1]
        S = compute_entropy(v, nbins, v_min, v_max, spwt)
        S_total_t[i] += S  # accumulate for total entropy
        delta_S.append(S - S0)

    plt.plot(ts[:DATA_TS_PHASE], delta_S, label=species_name)

# --------------------- Total Entropy Change Plot ---------------------
delta_S_total = [S - S_total_t0 for S in S_total_t]
plt.plot(ts[:DATA_TS_PHASE], delta_S_total, '--k', label="Total", linewidth=2)

# --------------------- Final Plot Formatting ---------------------
plt.xlabel(r'Time ($\omega_{pi} t$)')
plt.ylabel(r'$\Delta S = S(t) - S(0)$')
plt.title('Entropy Change vs Time for All Species')
plt.legend()
plt.grid(True)
plt.tight_layout()

# --------------------- Save and Show ---------------------
plt.savefig(pjoin(path_fig, "entropy_change_all_species.png"))
plt.show()
