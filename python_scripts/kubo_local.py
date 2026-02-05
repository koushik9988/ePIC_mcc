"""
Simple python script to analyze Electric field to study turbulence.
@kaushik kalita
@Date : 11/11/2025
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.constants import value as constants
import os
from tqdm import tqdm 

# -------------------- CONSTANTS --------------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# -------------------- INPUT CHECK --------------------
if len(sys.argv) != 3:
    print("Usage: python Dth_auto_full.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
file_name = 'result.h5'

output_dir = 'acf_plots'+'_'+particle_type
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

# -------------------- HDF5 OPEN --------------------
f = h5py.File(pjoin(path, file_name), 'r')

# -------------------- READ METADATA --------------------
metadata = f['/metadata']

NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
write_interval_phase = metadata.attrs['write_int_phase']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
normscheme = metadata.attrs['norm_scheme']
density = metadata.attrs['density']

species_metadata_path = f'/metadata_species/{particle_type}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']

q_mass_ratio = e / species_mass
print(species_spwt)

# -------------------- NORMALIZATION --------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi / wpe
if normscheme == 1 or normscheme == 2:
    mfactor = 1

DT = DT_coeff * (1.0 / wpe)

L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi

efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12

# -------------------- LOAD E-FIELD --------------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape

x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]

# -------------------- TIME ARRAY --------------------
time_full = np.arange(Nt) * DT_coeff * write_interval * mfactor
dt_field = DT_coeff * write_interval * mfactor

# -------------------- FFT (k) --------------------
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)
k_modes = kfft[:Nx // 2]
F_k_pos = F_k[:, :Nx // 2]

power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)
max_idx = peaks[np.argmax(power_spectrum[peaks])]

k_dominant = k_modes[max_idx]
E_k_t = F_k_pos[:, max_idx]

# -------------------- PHASE VELOCITY --------------------
dt_snapshot = time_full[1] - time_full[0]

omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)

omega_r = omega_freqs[:Nt // 2][np.argmax(np.abs(omega_spectrum[:Nt // 2]))]
v_phase = omega_r / k_dominant

print(f"Dominant k: {k_dominant:.4f}, Frequency: {omega_r:.4f}, V_phase: {v_phase:.4f}")

# -------------------- BOUNCE FREQUENCY WINDOW --------------------
tstart_idx = 1000
tau = 3900
tend_idx = tstart_idx + tau

start_time = max(0, int(round(tstart_idx / dt_field)))
end_time = min(Nt, int(round(tend_idx / dt_field)))

EF_slice = EF[tstart_idx:tend_idx, :]
Nt1, Nx1 = EF_slice.shape

print("shape of EF after slicing:", EF.shape)

x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]

# -------------------- FFT OVER SLICE --------------------
Fk = np.fft.fft(EF_slice, axis=1, norm='ortho')
kfft = 2.0 * np.pi * np.fft.fftfreq(Nx1, d=dx)

k_pos = kfft[:Nx1 // 2]
Fk_pos = Fk[:, :Nx1 // 2]

time_full = np.arange(Nt1) * dt_snapshot

Fk_pos = Fk_pos * efield_norm * (2.0 / np.sqrt(Nx))

Ek2 = np.abs(Fk_pos)**2

# -------------------- PEAK FINDING --------------------
power_spectrum = np.mean(Ek2, axis=0)
threshold = np.max(power_spectrum) * 0.001

peaks, _ = find_peaks(power_spectrum, height=threshold, distance=2)

N_modes = min(100, len(peaks))
top_order = np.argsort(power_spectrum[peaks])[-N_modes:]
top_idx = peaks[top_order]

k_dom = k_pos[top_idx] / L

k_avg = np.average(k_dom)
k_max = k_pos[top_idx[np.argmax(power_spectrum[top_idx])]] / L

E2_mean_mode = np.mean(Ek2[:, top_idx], axis=0)
k_modes_eff = k_pos[top_idx] / L

k_eff = np.sum(k_modes_eff * E2_mean_mode) / np.sum(E2_mean_mode)

print("k_eff =", k_eff)
print("k_dominant =", k_max)
print("k_avg:", k_avg)

# -------------------- BOUNCE FREQUENCY --------------------
Ek2_dom = Ek2[:, top_idx]
sum_Ek2_dom = np.sum(Ek2_dom, axis=1)
sum_Ek2_dom_avg = np.mean(sum_Ek2_dom)

E_eff = np.sqrt(sum_Ek2_dom_avg)
omega_b = np.sqrt((e * k_eff * E_eff) / species_mass)
omega_b_norm = omega_b / W

print("omega_b / wpi =", omega_b_norm)
print("Tb_effective =", 2 * np.pi / omega_b_norm)

Tb_eff = 2.0 * np.pi / omega_b_norm

# -------------------- TIME WINDOW FOR PARTICLES --------------------
E_k_amp = np.abs(E_k_t) * efield_norm * (2 / np.sqrt(Nx))
field_idx_peak = np.argmax(E_k_amp)

threshold = 0.1 * E_k_amp[field_idx_peak]
valid_indices = np.where(E_k_amp[field_idx_peak:] < threshold)[0]

if len(valid_indices) > 0:
    field_idx_end = field_idx_peak + valid_indices[0]
else:
    field_idx_end = len(E_k_amp) - 1

t1_time = 0.0
t2_time = tstart_idx

dt_particle = DT_coeff * write_interval_phase * mfactor

i1 = max(0, int(round(tstart_idx / dt_particle)))
i2 = min(Nt, int(round(tend_idx / dt_particle)))

print(f"Analysis Window Indices (Phase Snapshots): i1={i1}, i2={i2}")

lags_time_window = np.arange(i2 - i1) * dt_particle
num_lags = len(lags_time_window)

# -------------------- VELOCITY BINNING --------------------
particle_group = f[f"particle_{particle_type}"]
phasespace_keys = sorted(particle_group.keys(), key=int)

Nt_particle = len(phasespace_keys)
field_indices = np.round(np.linspace(0, Nt - 1, Nt_particle)).astype(int)

print("Scanning particle velocities for binning...")

num_particles_to_scan = min(species_num, 100000)
sample_velocities = []

for p in range(num_particles_to_scan):
    v = particle_group[phasespace_keys[i1]][p, 1]
    sample_velocities.append(v)

v_flat = np.array(sample_velocities)
vmin, vmax = np.percentile(v_flat, [5, 95])

num_bins = 50
bin_edges = np.linspace(vmin, vmax, num_bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

print(f"Velocity Bins: {num_bins} bins from {vmin:.4f} to {vmax:.4f}")

bin_acf_sum = np.zeros((num_bins, num_lags))
bin_counts = np.zeros(num_bins)

particle_start, particle_end = 0, num_particles_to_scan
if particle_end > species_num:
    particle_end = species_num

print("Starting Ensemble Accumulation...")

# -------------------- PROCESS PARTICLES --------------------
for p in tqdm(range(particle_start, particle_end), desc="Processing Particles"):

    v_start = particle_group[phasespace_keys[i1]][p, 1]

    if v_start < vmin or v_start >= vmax:
        continue

    bin_idx = int(np.digitize(v_start, bin_edges)) - 1

    if bin_idx < 0 or bin_idx >= num_bins:
        continue

    traj_window_keys = phasespace_keys[i1:i2]
    particle_x = np.array([particle_group[k][p, 0] for k in traj_window_keys])

    Nt_w = len(particle_x)
    E_traj = np.zeros(Nt_w)

    window_field_indices = field_indices[i1:i2]

    for j in range(Nt_w):
        f_idx = window_field_indices[j]

        xp = particle_x[j]
        i_x = int(np.floor(xp))
        dxp = xp - i_x

        i_x_next = (i_x + 1) % Nx

        E_traj[j] = EF[f_idx, i_x] * (1 - dxp) + EF[f_idx, i_x_next] * dxp

    E_env = (E_traj - np.mean(E_traj)) * efield_norm

    ac = np.correlate(E_env, E_env, mode='full')[Nt_w - 1:]

    if ac[0] != 0:
        ac = ac / ac[0]
    else:
        ac = np.zeros_like(ac)

    bin_acf_sum[bin_idx] += ac
    bin_counts[bin_idx] += 1

# -------------------- ENSEMBLE KUBO RESULTS --------------------
print("\nCalculating Ensemble Kubo Numbers...")

def exp_decay(t, tau, c):
    return np.exp(-t / tau) + c

bin_ku_list = []
bin_tau_list = []
valid_bins_vel = []

txt_out = pjoin(path_fig, f"ensemble_kubo_vs_vel_{particle_type}.txt")
f_out = open(txt_out, 'w')
f_out.write("# v_bin_center  <Ku>  tau_c  count\n")

for b in range(num_bins):

    count = bin_counts[b]
    if count < 10:
        continue

    avg_acf = bin_acf_sum[b] / count

    fit_end = int(0.20 * num_lags)

    try:
        popt, _ = curve_fit(exp_decay,lags_time_window[:fit_end],avg_acf[:fit_end],p0=(dt_particle * 5, 0.0),maxfev=5000)
        tau_c = abs(popt[0])
    except:
        tau_c = 0.0

    ac = avg_acf  # rename for clarity

    # Determine where ACF crosses zero
    if np.any(ac <= 0):
        cutoff_i = np.where(ac <= 0)[0][0]   # first zero-crossing index
    else:
        cutoff_i = len(ac)                   # no zero crossing → integrate full range

    # Integral time scale (Wikipedia definition)
    tau_int = np.trapz(ac[:cutoff_i], lags_time_window[:cutoff_i])

    #Ku = tau_int / Tb_eff
    Ku = tau_c / Tb_eff 


    bin_ku_list.append(Ku)
    bin_tau_list.append(tau_c)
    valid_bins_vel.append(bin_centers[b])

    f_out.write(f"{bin_centers[b]:.6e} {Ku:.6e} {tau_c:.6e} {int(count)}\n")


    #####plots
    # ---------------------------------------------------------
    plt.figure(figsize=(7, 5))

    # Plot ACF
    plt.plot(lags_time_window, avg_acf, 'b-', label='Ensemble ACF')

    # Plot exponential fit curve
    fit_curve = np.exp(-lags_time_window / tau_c) + popt[1]
    plt.plot(lags_time_window, fit_curve, 'r--', label='Exp Fit')

    plt.title(f"Bin {b}: v = {bin_centers[b]:.3f}, count = {int(count)}")
    plt.xlabel("Time Lag")
    plt.ylabel("ACF")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Save figure
    acf_plot_name = pjoin(path_fig, f"acf_bin_{b:03d}_v{bin_centers[b]:.3f}.png")
    plt.savefig(acf_plot_name, dpi=150)
    plt.close()
    ####plots

f_out.close()

print(f"Saved Ensemble Data to {txt_out}")

# -------------------- PLOT --------------------
plt.figure(figsize=(10, 6))
plt.plot(valid_bins_vel, bin_ku_list, 'o-', linewidth=2, label='Ensemble Averaged Kubo')

plt.xlabel(r"$v$")
plt.ylabel(r"$Ku = \tau_{ac} / T_b$")
plt.legend()
plt.grid(True, alpha=0.3)

plot_png = pjoin(path_fig, f"ensemble_kubo_{particle_type}.png")
plt.savefig(plot_png, dpi=200)

print(f"Saved Plot to {plot_png}")
plt.show()
