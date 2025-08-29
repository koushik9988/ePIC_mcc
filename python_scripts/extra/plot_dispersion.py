import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.constants import value as constants
import h5py
import sys
from os.path import join as pjoin

# Constants
eps0 = constants("electric constant")
kb = constants("Boltzmann constant")
e = constants("elementary charge")
me = constants("electron mass")
AMU = constants("atomic mass constant")

# Load input
if len(sys.argv) != 2:
    print("Usage: python3 plot_dispersion_generalized.py <data_folder>")
    sys.exit(1)

data_path = sys.argv[1]
file_path = pjoin(data_path, 'result.h5')

try:
    f = h5py.File(file_path, 'r')
except FileNotFoundError:
    print(f"Error: HDF5 file not found at {file_path}")
    sys.exit(1)

# Read global metadata
meta = f["/metadata"]
norm_scheme = meta.attrs.get("norm_scheme", "unnormalized")
wpe = meta.attrs.get("wpe", None)
LDe = meta.attrs.get("LDe", None)

if wpe is None or LDe is None:
    print("Error: Missing 'wpe' or 'LDe' in metadata")
    sys.exit(1)

# Debug: Print raw metadata
print(f"Metadata: wpe={wpe:.3e}, LDe={LDe:.3e}, norm_scheme={norm_scheme}")

# Normalization factors
L = LDe              # Length normalization (Debye length)
W = wpe              # Frequency normalization (plasma frequency)
V = L * W            # Velocity normalization

# Load species info
species_order = f["/metadata_species"].attrs.get("species_order", None)
if species_order is None:
    print("Error: Missing 'species_order' attribute in /metadata_species")
    sys.exit(1)

species_list = []
for name in species_order:
    try:
        sp = f[f"/metadata_species/{name}"].attrs
        Tj = sp["temperature"]
        mj = sp["mass"]
        nj = sp["density"]
        v0j = sp.get("streaming_velocity", 0.0)

        # Debug: Print raw species data
        print(f"Raw data for {name}: density={nj:.3e}, mass={mj:.3e}, temperature={Tj:.3e}, v0={v0j:.3e}")

        # Normalized quantities
        vth = np.sqrt(kb * Tj / mj) / V  # Thermal velocity
        wp = np.sqrt(nj * e**2 / (eps0 * mj)) / W  # Plasma frequency
        v0 = v0j / V  # Streaming velocity

        # Skip species with zero density or invalid wp
        if nj <= 0 or not np.isfinite(wp):
            print(f"Warning: Skipping species {name} due to invalid density ({nj:.3e}) or wp ({wp:.3e})")
            continue

        species_list.append({
            "name": name,
            "vth": vth,
            "wp": wp,
            "v0": v0,
        })
    except KeyError as e:
        print(f"Error: Missing attribute {e} for species {name}")
        sys.exit(1)

f.close()

# Debug: Print species info
print(f"Loaded {len(species_list)} species:")
for sp in species_list:
    print(f"  {sp['name']}: wp={sp['wp']:.3e}, vth={sp['vth']:.3e}, v0={sp['v0']:.3e}")

if not species_list:
    print("Error: No valid species to process. Exiting.")
    sys.exit(1)

# Dispersion relation function
def dispersion_func(omega_complex, k, damping=1e-10):
    omega = omega_complex[0] + 1j * omega_complex[1]
    total = 0.0
    for sp in species_list:
        wp = sp["wp"]
        vth = sp["vth"]
        v0 = sp["v0"]
        # Add small damping to denominator to avoid division by zero
        denom = (omega - k * v0)**2 - (k * vth)**2 + damping
        total += wp**2 / denom
    diff = 1 - total
    return [diff.real, diff.imag]

# Solve dispersion relation
k_vals = np.linspace(0.01, 2.0, 200)
omega_real = []
omega_imag = []
guess = [1.0, 0.001]  # More robust initial guess
success_count = 0

for k in k_vals:
    try:
        sol = root(dispersion_func, guess, args=(k,), method='hybr', options={'xtol': 1e-8, 'maxfev': 1000})
        if sol.success:
            w = sol.x[0] + 1j * sol.x[1]
            omega_real.append(w.real)
            omega_imag.append(w.imag)
            guess = sol.x  # Update guess for next iteration
            success_count += 1
        else:
            omega_real.append(np.nan)
            omega_imag.append(np.nan)
            print(f"Warning: Solver failed for k={k:.3f}")
    except Exception as e:
        print(f"Error: Solver exception for k={k:.3f}: {e}")
        omega_real.append(np.nan)
        omega_imag.append(np.nan)

# Debug: Check solver success rate
print(f"Solver success rate: {success_count}/{len(k_vals)} ({success_count/len(k_vals)*100:.1f}%)")
print(f"Valid omega_real points: {np.sum(~np.isnan(omega_real))}/{len(omega_real)}")
print(f"Valid omega_imag points: {np.sum(~np.isnan(omega_imag))}/{len(omega_imag)}")

# Plotting with ax objects
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Check if data is valid
if np.sum(~np.isnan(omega_real)) > 0 and np.sum(~np.isnan(omega_imag)) > 0:
    # Real part plot
    ax1.plot(k_vals, omega_real, label="Re($\omega$)", color='blue')
    ax1.set_xlabel(r"$k \lambda_D$")
    ax1.set_ylabel(r"Re($\omega/\omega_{pe}$)")
    ax1.set_title(f"Real Part (norm: {norm_scheme})")
    ax1.grid(True)
    ax1.legend()

    # Imaginary part plot
    ax2.plot(k_vals, omega_imag, label="Im($\omega$)", color='red')
    ax2.set_xlabel(r"$k \lambda_D$")
    ax2.set_ylabel(r"Im($\omega/\omega_{pe}$)")
    ax2.set_title(f"Imaginary Part (norm: {norm_scheme})")
    ax2.grid(True)
    ax2.legend()
else:
    print("Warning: No valid data to plot. Displaying cold plasma dispersion relation.")
    # Fallback: Plot cold plasma dispersion relation (wp^2 = sum(wp_j^2))
    omega_cold = np.sqrt(sum(sp["wp"]**2 for sp in species_list) + k_vals**2)
    ax1.plot(k_vals, omega_cold, label="Cold Plasma Re($\omega$)", color='blue', linestyle='--')
    ax1.set_xlabel(r"$k \lambda_D$")
    ax1.set_ylabel(r"Re($\omega/\omega_{pe}$)")
    ax1.set_title(f"Cold Plasma Approx (norm: {norm_scheme})")
    ax1.grid(True)
    ax1.legend()

    ax2.plot(k_vals, np.zeros_like(k_vals), label="Cold Plasma Im($\omega$)", color='red', linestyle='--')
    ax2.set_xlabel(r"$k \lambda_D$")
    ax2.set_ylabel(r"Im($\omega/\omega_{pe}$)")
    ax2.set_title(f"Cold Plasma Approx (norm: {norm_scheme})")
    ax2.grid(True)
    ax2.legend()

fig.tight_layout()
fig.suptitle(f"Generalized Dispersion Relation ({len(species_list)} species)", y=1.05)
plt.show(block=True)