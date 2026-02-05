import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python3 diffusion.py <path_to_data> <particle_type>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]

plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Load HDF5 data
file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

# Metadata
metadata = f['/metadata'].attrs
NUM_TS = metadata['NUM_TS']
write_interval_phase = metadata['write_int_phase']
DT_coeff = metadata['DT_coeff']
wpi = metadata['wpi']
wpe = metadata['wpe']
LDe = metadata['LDe']
LDi = metadata['LDi']

vthe = LDe * wpe
vthi = LDi * wpi

mFactor = wpi / wpe

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

print(f"Total number of phase space data: {DATA_TS_PHASE}")

data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mFactor  # ω_pi * t

# Reference velocity: initial mean (t=0)
data_phase0 = f[f"particle_{particle_type}/0"]
v0 = data_phase0[:, 1]
v0_mean = np.mean(v0)
print("v_mean =", v0_mean)

#empty lists to store results
sigma_v2 = []
t_vals = []

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    data_phase_t = f[f"particle_{particle_type}/{j}"]
    v_t = data_phase_t[:, 1]
    sigma2 = np.mean((v_t - v0_mean) ** 2)
    sigma_v2.append(sigma2)
    t_vals.append(j*DT_coeff*mFactor)

sigma_v2 = np.array(sigma_v2)
t_vals = np.array(t_vals)


plt.figure()
plt.plot(t_vals, sigma_v2, label=r'$\sigma_v^2(t)$')
plt.legend()
plt.xlabel(r'$\omega_{pi} t$')
plt.ylabel(r'$(v(t) - v(0)^2$')
plt.savefig(pjoin(path_fig, f"velocity_diffusion_{particle_type}.png"))
plt.show()



#staring point for fit
t_start = 150

#Filter out all points before t_start
mask = t_vals >= t_start
t_fit = t_vals[mask]
sigma_fit = sigma_v2[mask]


coeff = np.polyfit(t_fit, sigma_fit, 1)
fit_line = np.polyval(coeff, t_fit)

# Plot
plt.plot(t_vals, sigma_v2, label=r'$\sigma_v^2(t)$')
plt.plot(t_fit, fit_line, 'r--', label=f'slope={coeff[0]:.3e}')

plt.legend()
plt.xlabel(r'$\omega_{pi} t$')
plt.ylabel(r'$(v(t) - v(0))^2$')
plt.savefig(pjoin(path_fig, f"velocity_diffusion_{particle_type}.png"))
plt.show()

print("diffusion coefficient D =", coeff[0]/2)

"""
# Linear fit 
# numpy.polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False)
coeff = np.polyfit(t_vals, sigma_v2, 1)
fit_line = np.polyval(coeff, t_vals)
plt.plot(t_vals, sigma_v2, label=r'$\sigma_v^2(t)$')
plt.plot(t_vals, fit_line, 'r--', label='Linear fit')

plt.legend()
plt.xlabel(r'$\omega_{pi} t$')
plt.ylabel(r'$(v(t) - v(0)^2$')
plt.savefig(pjoin(path_fig, f"velocity_diffusion_{particle_type}.png"))
plt.show()
"""