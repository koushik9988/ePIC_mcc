import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
import os
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks

if len(sys.argv) != 3:
    print("Usage: python3 mode_diffusion.py <path_to_data> <particle_type>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]

plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

#Metadata 
metadata = f['/metadata'].attrs
NUM_TS = metadata['NUM_TS']
write_interval_phase = metadata['write_int_phase']
DT_coeff = metadata['DT_coeff']
wpi = metadata['wpi']
wpe = metadata['wpe']
LDe = metadata['LDe']
LDi = metadata['LDi']
mFactor = wpi / wpe  # For time scaling to ω_pi t
write_interval = metadata['write_int']  # For fields
vthe = LDe * wpe
vthi = LDi * wpi
q_over_m = 1.0 # Charge to mass ratio, set to 1 for normalized units

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1
print(f"Total particle snapshots: {DATA_TS_PHASE}")


x_part = []
vx_part = []
t_part = []  # ω_pi t
data_phase0 = f[f"particle_{particle_type}/0"]
v0_mean = np.mean(data_phase0[:, 1])
Np = len(data_phase0)

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data_phase_t = f[f"particle_{particle_type}/{j}"]
    x_part.append(data_phase_t[:, 0])  # x
    vx_part.append(data_phase_t[:, 1])  # vx
    t_part.append(j * DT_coeff * mFactor)  # ω_pi t

x_part = np.array(x_part).T  # (Np, Nt)
vx_part = np.array(vx_part).T  # (Np, Nt)
t_part = np.array(t_part)     # (Nt,)
Nt_part = len(t_part)
dt_part = t_part[1] - t_part[0] if Nt_part > 1 else 1.0

print(f"Particle data shape: x,vx (Np={Np}, Nt={Nt_part})")

# Load field data: E(x,t)
electric_field_data = []
time_steps_field = sorted(map(int, f['fielddata/efield'].keys()))
for ts_f in time_steps_field:
    EF_data = f[f'fielddata/efield/{ts_f}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # (Nt_field, Nx)
NC = metadata['NC']
x_grid = np.linspace(0, NC, EF.shape[1])
dx = x_grid[1] - x_grid[0]
DT_field = DT_coeff * (1.0 / wpe)
t_field = np.arange(EF.shape[0]) * DT_coeff * write_interval * mFactor  # Scale to ω_pi t
Nt_field = len(t_field)
dt_field = t_field[1] - t_field[0] if Nt_field > 1 else 1.0

print(f"Field data shape: E (Nt_field={Nt_field}, Nx={NC})")

# Spatial FFT for fields
kfft = 2 * np.pi * fftfreq(NC, dx)
k_modes = kfft[:NC // 2]
F_k_pos = fft(EF, axis=1, norm='ortho')[:, :NC // 2]
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)

# Find peaks (modes)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.001, distance=2)
k_peaks = k_modes[peaks]
print("Detected modes k:", k_peaks)

epsilon = 0.1  # Broadening for analytic delta
Delta_t_idx = min(5, Nt_part // 10)  # Short lag for traj
t_start_fit = 10  # For total D comparison

# Method 1: Analytic Quasilinear D_k (avg over particles, assume ω_k=0)
def broadened_delta(kv, eps=epsilon):
    return np.exp(-kv**2 / (2 * eps**2)) / (eps * np.sqrt(2 * np.pi))

D_k_analytic = np.zeros(len(k_modes))
for ik, kk in enumerate(k_modes):
    if abs(kk) < 1e-6: continue
    P_k = power_spectrum[ik]
    kv_res = kk * vx_part.flatten()
    resonance_avg = np.mean(broadened_delta(kv_res))
    D_k_analytic[ik] = (np.pi * q_over_m**2 / abs(kk)) * P_k * resonance_avg

# Restrict to peaks for comparison
D_k_analytic_peaks = D_k_analytic[peaks]
total_D_analytic = np.sum(D_k_analytic)

# Method 2: Trajectory-Based D_k
D_k_traj = np.zeros(len(k_modes))
for ik, kk in enumerate(k_modes):
    if abs(kk) < 1e-6 or Nt_part < 2: continue
    phase = np.exp(-1j * kk * x_part)  # (Np, Nt)
    v_tilde_k = vx_part * phase
    v_tilde_k_shift = np.roll(v_tilde_k, -Delta_t_idx, axis=1)
    delta_v_tilde_k = v_tilde_k_shift - v_tilde_k
    msd = np.mean(np.abs(delta_v_tilde_k[:, :-Delta_t_idx])**2)
    D_k_traj[ik] = msd / (2 * Delta_t_idx * dt_part)

D_k_traj_peaks = D_k_traj[peaks]
total_D_traj = np.sum(D_k_traj)

# Total diffusion from particles
sigma_v2 = []
for i in range(Nt_part):
    v_t = vx_part[:, i]
    sigma2 = np.mean((v_t - v0_mean)**2)
    sigma_v2.append(sigma2)
sigma_v2 = np.array(sigma_v2)
mask = t_part >= t_start_fit
t_fit = t_part[mask]
sigma_fit = sigma_v2[mask]
if len(t_fit) > 1:
    coeff = np.polyfit(t_fit, sigma_fit, 1)
    total_D_empirical = coeff[0] / 2
else:
    total_D_empirical = 0.0
print(f"Empirical total D: {total_D_empirical:.3e}")
print(f"Analytic total D: {total_D_analytic:.3e}")
print(f"Traj total D: {total_D_traj:.3e}")

# Plots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# 1. Power spectrum
axs[0,0].plot(k_modes, power_spectrum, label='Power')
axs[0,0].semilogy()
axs[0,0].plot(k_peaks, power_spectrum[peaks], 'ro', label='Peaks')
axs[0,0].set_xlabel('k'); axs[0,0].set_ylabel('|E_k|^2'); axs[0,0].legend(); axs[0,0].grid()

# 2. D_k comparison per mode (peaks)
axs[0,1].semilogy(k_peaks, D_k_analytic_peaks, 'b-', label='Analytic')
axs[0,1].semilogy(k_peaks, D_k_traj_peaks, 'r--', label='Traj')
axs[0,1].set_xlabel('k'); axs[0,1].set_ylabel('D_k'); axs[0,1].legend(); axs[0,1].grid()

# 3. sigma_v^2(t)
axs[1,0].plot(t_part, sigma_v2, label=r'$\sigma_v^2(t)$')
if len(t_fit) > 1:
    fit_line = np.polyval(coeff, t_fit)
    axs[1,0].plot(t_fit, fit_line, 'r--', label=f'Empirical D={total_D_empirical:.3e}')
axs[1,0].set_xlabel(r'$\omega_{pi} t$'); axs[1,0].set_ylabel(r'$\sigma_v^2$'); axs[1,0].legend(); axs[1,0].grid()

# 4. Total D comparison bar
methods = ['Empirical', 'Analytic', 'Traj']
D_totals = [total_D_empirical, total_D_analytic, total_D_traj]
axs[1,1].bar(methods, D_totals)
axs[1,1].set_ylabel('Total D'); axs[1,1].grid()

plt.tight_layout()
plt.savefig(pjoin(path_fig, f"mode_diffusion_comparison_{particle_type}.png"), dpi=300)
plt.show()

# Print per-mode comparison
print("\nPer-mode D_k comparison (for peaks):")
for i, k in enumerate(k_peaks):
    print(f"k={k:.2f}: Analytic={D_k_analytic_peaks[i]:.3e}, Traj={D_k_traj_peaks[i]:.3e}, Ratio={D_k_analytic_peaks[i]/D_k_traj_peaks[i]:.2f}")