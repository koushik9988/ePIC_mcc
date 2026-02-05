import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks

# ---------------------------------------------------------------- #
# Usage:
# python plot_phase_space_kubo_gt1.py <path_to_simulation> <path_to_csv_folder> <species_name>
# ---------------------------------------------------------------- #

if len(sys.argv) != 4:
    print("Usage: python plot_phase_space_kubo_gt1.py <sim_data_path> <csv_folder> <species_name>")
    sys.exit(1)

data_path = sys.argv[1]
csv_folder = sys.argv[2]
species_name = sys.argv[3]

# ---------------------------------------------------------------- #
# LOAD KUBO CSV
# ---------------------------------------------------------------- #
csv_file = pjoin(csv_folder, f"per_particle_kubo_{species_name}.csv")
df = pd.read_csv(csv_file)

df.columns = df.columns.str.strip()
df["kubo"] = pd.to_numeric(df["kubo"], errors="coerce")

particles_gt1 = df[df["kubo"] > 0.3]["particle_index"].astype(int).tolist()

#particles_gt1 = df[(df["kubo"] >= 0.9) & (df["kubo"] <= 1)]["particle_index"].astype(int).tolist()


print(f"\nFound {len(particles_gt1)} particles with Kubo > 1")
print("Particle indices:", particles_gt1)

if len(particles_gt1) == 0:
    print("No particles with Kubo > 1. Exiting.")
    sys.exit(0)


file_name = 'result.h5'
f = h5py.File(pjoin(data_path, file_name), 'r')

meta = f['/metadata']
NC = meta.attrs['NC']
NUM_TS = meta.attrs['NUM_TS']
write_interval = meta.attrs['write_int']
write_interval_phase = meta.attrs['write_int_phase']
DT_coeff = meta.attrs['DT_coeff']
wpe = meta.attrs['wpe']
wpi = meta.attrs['wpi']
norm_scheme = meta.attrs['norm_scheme']

mfactor = (wpi / wpe) if norm_scheme not in [1, 2] else 1


electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for t in time_steps:
    E = f[f'fielddata/efield/{t}'][:]
    electric_field_data.append(E)

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x_grid = np.linspace(0, NC, Nx)
dx = x_grid[1] - x_grid[0]

dt_field = DT_coeff * write_interval * mfactor
F_k = np.fft.fft(EF, axis=1, norm='ortho')

kfft = 2*np.pi*np.fft.fftfreq(Nx, dx)
k_modes = kfft[:Nx//2]
F_k_pos = F_k[:, :Nx//2]

power = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power, height=np.max(power)*0.04, distance=2)

k_dominant = k_modes[peaks[np.argmax(power[peaks])]]
E_k_t = F_k_pos[:, peaks[np.argmax(power[peaks])]]

omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2*np.pi*np.fft.fftfreq(Nt, dt_field)

omega_pos = omega_freqs[:Nt//2]
omega_spectrum_pos = omega_spectrum[:Nt//2]

omega_r = omega_pos[np.argmax(np.abs(omega_spectrum_pos))]
v_phase = omega_r / k_dominant

print(f"\nDominant phase velocity v_phase = {v_phase:.4f}")



particle_group = f[f'particle_{species_name}']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1
time_arr = np.arange(DATA_TS_PHASE) * write_interval_phase * DT_coeff * mfactor

plt.figure(figsize=(8,6))

for pidx in particles_gt1:

    x_traj = []
    v_traj = []

    for i in range(DATA_TS_PHASE):
        ts = i * write_interval_phase
        pdata = particle_group[str(ts)]
        x_traj.append(pdata[pidx, 0])
        v_traj.append(pdata[pidx, 1])

    x_traj = np.array(x_traj)
    v_traj = np.array(v_traj)

    # unwrap
    x_unwrapped = np.unwrap(x_traj * 2*np.pi / NC) * NC / (2*np.pi)

    # wave frame
    x_wave = x_unwrapped - v_phase * time_arr
    v_wave = v_traj

    plt.plot(x_wave, v_wave, lw=1.2)#, label=f"p={pidx}")

plt.xlabel(r"$x - v_p t$")
plt.ylabel(r"$v$")
plt.title(f"Wave-Frame Phase-Space\nAll Kubo > 1 Particles ({species_name})")
plt.grid(True)


#plt.legend(fontsize=8, ncol=2)

plt.tight_layout()
plt.show()

f.close()
