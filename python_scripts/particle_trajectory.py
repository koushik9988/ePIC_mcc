import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.constants import value as constants

# --------------------- Argument Parsing ---------------------
if len(sys.argv) != 4:
    print("Usage: python3 track_particle.py <path_to_data> <particle_type> <particle_index>")
    sys.exit(1)

data_path = sys.argv[1]
particle_type = sys.argv[2]
particle_index = int(sys.argv[3])

file_name = 'result.h5'
file_path = pjoin(data_path, file_name)

# --------------------- Load File & Metadata ---------------------
f = h5py.File(file_path, 'r')
meta = f['/metadata']

NC = meta.attrs['NC']
NUM_TS = meta.attrs['NUM_TS']
write_interval_phase = meta.attrs['write_int_phase']
wpe = meta.attrs['wpe']
wpi = meta.attrs['wpi']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1
mFactor = wpi / wpe

# Time array in ω_pi * t
ts = f['time_var/kinetic_energy'][:, 0]
ts = ts * mFactor

# --------------------- Track Particle ---------------------
x_traj = []
v_traj = []
time_arr = []

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data = f[f"particle_{particle_type}/{j}"]
    x_traj.append(data[particle_index, 0])
    v_traj.append(data[particle_index, 1])
    time_arr.append(ts[i])

f.close()

x_traj = np.array(x_traj)
v_traj = np.array(v_traj)
time_arr = np.array(time_arr)

# --------------------- Plot ---------------------
fig, (ax1, ax2 ,ax3) = plt.subplots(3, 1, figsize=(8, 6))

ax1.plot(time_arr, x_traj, label='x(t)')
ax1.set_ylabel('$x$')
ax1.grid(True)
ax1.legend()

ax2.plot(time_arr, v_traj, color='orange', label='v(t)')
ax2.set_xlabel('[$\omega_{pi} t$]')
ax2.set_ylabel('$v$')
ax2.grid(True)
ax2.legend()


ax3.plot(x_traj, v_traj, color='orange', label='$x-v$')
ax3.set_xlabel('$x$')
ax3.set_ylabel('$v$')
ax3.grid(True)
ax3.legend()

fig.suptitle(f'Trajectory of Particle #{particle_index} ({particle_type})')
#plt.tight_layout()
plt.show()

