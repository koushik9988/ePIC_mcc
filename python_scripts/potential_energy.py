import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin, basename

# Check usage
if len(sys.argv) < 2:
    print("Usage: python potnetial_energy.py <path1> <path2> ...")
    sys.exit(1)

paths = sys.argv[1:]

plt.figure()

for path in paths:
    file_name = 'result.h5'
    file_path = pjoin(path, file_name)

    if not os.path.isfile(file_path):
        print(f"Skipping: {file_path} (not found)")
        continue

    # Output directory
    output_dir = pjoin(path, "plots")
    os.makedirs(output_dir, exist_ok=True)

    # Open HDF5
    f = h5py.File(file_path, 'r')

    metadata_group = f['/metadata']
    wpe = metadata_group.attrs['wpe']
    wpi = metadata_group.attrs['wpi']
    normscheme = metadata_group.attrs['norm_scheme']

    mfactor = wpi / wpe

    data = f["time_var/kinetic_energy"]
    ts = data[:, 0] * mfactor

    spno = metadata_group.attrs["spno"]
    pe = data[:, 3 * spno + 1]

    label_name = basename(os.path.normpath(path))

    plt.plot(ts, pe, label=f"{label_name}", linewidth=1.2)

    f.close()

# Formatting
plt.yscale("log")
plt.ylim([1e-4, 1e0])
plt.xlabel(r"$\omega_{pe}t$")
plt.ylabel("Potential Energy")
plt.ticklabel_format(style='sci', axis='y')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()

plt.savefig("pe_multi.png", dpi=1200)
plt.show()
