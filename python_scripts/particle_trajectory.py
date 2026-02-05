import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks
from scipy.constants import value as constants

# ---------------------------------------------------------------- #
# Usage
if len(sys.argv) != 4:
    print("Usage: python scripts.py <path_to_data> <species_name> <particle_index>")
    sys.exit(1)

data_path = sys.argv[1]
species_name = sys.argv[2]
particle_index = int(sys.argv[3])

# ---------------------------------------------------------------- #
# Load main HDF5 file
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
LDe = meta.attrs['LDe']
LDi = meta.attrs['LDi']
norm_scheme = meta.attrs['norm_scheme']

mfactor = wpi / wpe
if norm_scheme in [1, 2]:
    mfactor = 1

# ---------------------------------------------------------------- #
# LOAD ELECTRIC FIELD AND COMPUTE PHASE VELOCITY (from trapped_width.py)

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
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)

k_modes = kfft[:Nx//2]
F_k_pos = F_k[:, :Nx//2]
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)

peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)
max_idx = peaks[np.argmax(power_spectrum[peaks])]
k_dominant = k_modes[max_idx]

E_k_t = F_k_pos[:, max_idx]

# temporal FFT
omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_field)

omega_pos_freqs = omega_freqs[:Nt//2]
omega_pos_spectrum = omega_spectrum[:Nt//2]

omega_peak_idx = np.argmax(np.abs(omega_pos_spectrum))
omega_r = omega_pos_freqs[omega_peak_idx]

v_phase = omega_r / k_dominant   # NORMALIZED PHASE VELOCITY

print(f"\n wave phase velocity v_phase = {v_phase:.4f})")

# ---------------------------------------------------------------- #

particle_group = f[f'particle_{species_name}']

time_arr = []
x_traj = []
v_traj = []

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

for i in range(DATA_TS_PHASE):
    ts_index = i * write_interval_phase
    p_data = particle_group[str(ts_index)]

    x_traj.append(p_data[particle_index, 0])
    v_traj.append(p_data[particle_index, 1])
    time_arr.append(i * write_interval_phase * DT_coeff * mfactor)

x_traj = np.array(x_traj)
v_traj = np.array(v_traj)
time_arr = np.array(time_arr)


x_test = np.unwrap(x_traj * 2*np.pi / NC) * NC / (2*np.pi)
x_wave = x_test - v_phase * time_arr


#x_wave = x_traj - v_phase * time_arr
v_wave = v_traj #- v_phase


plt.figure()
plt.plot(x_wave, v_wave)
plt.xlabel(r"$x-v_p t$")
plt.ylabel(r"$v$")
plt.title(f"Particle ID = {particle_index} Phase-space in Wave Frame kubo = 0.01")
plt.grid(True)
plt.tight_layout()
plt.show()

f.close()
