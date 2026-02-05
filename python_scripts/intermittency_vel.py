"""
@kaushik kalita
@Date : 02/12/2025
Index-based space-lagged velocity structure functions
"""

import numpy as np
import h5py
import sys
import os
import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.constants import value as constants

# ---------------- constants ----------------
eps0 = constants('electric constant')
kb   = constants('Boltzmann constant')
me   = constants('electron mass')
AMU  = constants('atomic mass constant')
e    = constants('elementary charge')
# -------------------------------------------

# ---------------- arguments ----------------
if len(sys.argv) != 3:
    print("Usage: python space_lag_index.py <path_to_data> <particle_type>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
file_name = 'result.h5'
# -------------------------------------------

# ---------------- output -------------------
output_dir = 'new_analysis'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)
# -------------------------------------------

# ---------------- open file ----------------
f = h5py.File(pjoin(path, file_name), 'r')

# ---------------- metadata -----------------
metadata = f['/metadata']
NUM_TS = metadata.attrs['NUM_TS']
write_interval_phase = metadata.attrs['write_int_phase']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
normscheme = metadata.attrs['norm_scheme']

species_metadata = f[f'/metadata_species/{particle_type}']
species_num = species_metadata.attrs['num_particles']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme in (1, 2):
    mfactor = 1.0

L = LDe
W = wpe
if normscheme in (2, 4):
    L = LDi
    W = wpi

DT = DT_coeff / W
dt_phase = write_interval_phase * DT_coeff * mfactor
# -------------------------------------------

# ---------------- load data ----------------
particle_group_name = f"particle_{particle_type}"

if particle_group_name not in f:
    print(f"Error: {particle_group_name} not found")
    sys.exit(1)

keys = sorted([int(k) for k in f[particle_group_name].keys()])

x_list  = []
vel_list = []

for j in keys:
    ds = f[f"{particle_group_name}/{j}"]
    x_list.append(ds[:, 0].astype(float))
    vel_list.append(ds[:, 1].astype(float))

x   = np.array(x_list)
vel = np.array(vel_list)

nt, npart = vel.shape
print(f"Loaded {nt} snapshots, {npart} particles")

# -------------------------------------------

# ---------- space-lag index function -------
def space_lag_index(x_t, v_t, di):
    """
    Index-based space-lagged velocity increment
    """
    idx = np.argsort(x_t)
    v_sorted = v_t[idx]
    dv = v_sorted[di:] - v_sorted[:-di]
    return dv
# -------------------------------------------

# ---------- structure function -------------
def structure_function(dv, p):
    return np.mean(np.abs(dv)**p)
# -------------------------------------------

# ---------- parameters ---------------------
lags = np.arange(1, npart // 20)   # safe range
orders = [2, 4]                    # S2, S4
t_skip = 5                          # time averaging stride
# -------------------------------------------

# ---------- compute structure functions ----
Sp = {p: np.zeros(len(lags)) for p in orders}
count = 0

for it in range(0, nt, t_skip):
    x_t = np.mod(x[it], L)          # periodic domain
    v_t = vel[it]

    for k, di in enumerate(lags):
        dv = space_lag_index(x_t, v_t, di)
        for p in orders:
            Sp[p][k] += structure_function(dv, p)

    count += 1

for p in orders:
    Sp[p] /= count

# -------------------------------------------

# ---------- save data ----------------------
np.savez(
    pjoin(path_fig, f"space_lag_index_{particle_type}.npz"),
    lags=lags,
    **{f"S{p}": Sp[p] for p in orders}
)
# -------------------------------------------

# ---------- plot ---------------------------
plt.figure(figsize=(6, 4))
for p in orders:
    plt.loglog(lags, Sp[p], marker='o', label=rf"$S_{p}$")

plt.xlabel(r"Index lag $\Delta i$")
plt.ylabel(r"$S_p(\Delta i)$")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.tight_layout()
plt.savefig(pjoin(path_fig, f"space_lag_structure_{particle_type}.png"), dpi=300)
plt.close()

print("Space-lagged structure functions computed and saved")
