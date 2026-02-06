import numpy as np
import h5py
import sys
import os
import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.optimize import curve_fit
from scipy.constants import value as constants
from tqdm import tqdm
from scipy.signal import find_peaks

# ---------------------- CLI ----------------------
if len(sys.argv) != 3:
    print("Usage: python kubo_ensemble_new.py <path_to_data> <particle_type>")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]

# ------------------ OUTPUT ----------------------
outdir = "kubo_ensemble_time"
path_out = pjoin(path, outdir)
os.makedirs(path_out, exist_ok=True)

file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

# ------------------- METADATA ----------------------
metadata_group = f['/metadata']

eps0 = constants('electric constant')
kb   = constants('Boltzmann constant')
me   = constants('electron mass')
e    = constants('elementary charge')

NC   = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Ti = metadata_group.attrs['Ti']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
normscheme = metadata_group.attrs['norm_scheme']
EV_TO_K = 11604.52

# ------------------- PLASMA PARAMS =====================
ni0 = n0
ne0 = n0 / (1 + alp + beta)

LD  = np.sqrt(eps0 * kb * Te * EV_TO_K / (ne0 * e**2))
LDi = np.sqrt(eps0 * kb * Ti * EV_TO_K / (ni0 * e**2))
wpe = np.sqrt(ne0 * e**2 / (eps0 * me))
wpi = np.sqrt(ni0 * e**2 / (eps0 * mi))

mfactor = wpi / wpe
if normscheme in (1, 2):
    mfactor = 1

efield_norm = (n0 * e * LDi) / eps0

dt_phase = write_interval_phase * DT_coeff * mfactor
dt_field = write_interval * DT_coeff * mfactor

# ------------------- LOAD FIELD --------------------==
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for ts in time_steps:
    electric_field_data.append(f[f'fielddata/efield/{ts}'][:])

EF = np.vstack(electric_field_data)
Nt_field, Nx = EF.shape

x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]


# -------------LOAD PARTICLES ----------------
particle_group = f[f"particle_{particle_type}"]
phasespace_keys = sorted(particle_group.keys(), key=int)
Nt_particle = len(phasespace_keys)

field_indices = np.round(np.linspace(0, Nt_field - 1, Nt_particle)).astype(int)

# ------------------- VELOCITY BINNING -------------------
num_particles_to_scan = min(nParticlesB, 200000)

sample_velocities = []

for p in range(num_particles_to_scan):
    v = particle_group[phasespace_keys[0]][p, 1]
    sample_velocities.append(v)

v_flat = np.array(sample_velocities)
vmin, vmax = np.percentile(v_flat, [5, 95])

nv_bins = 50
v_edges = np.linspace(vmin, vmax, nv_bins + 1)
v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])

# ------------------- TIME WINDOWS --------------------
#physical time in (wpi*t) units
tstart_phys =  int(dt_phase *  506)#80#600#690.0
tend_phys   = int(dt_phase *  1973)#320#2200.0

print("field indices are in wpi t unit",tstart_phys,"and",tend_phys)
print("checking",(  NUM_TS*DT_coeff*mfactor/dt_phase))

window_phys = 100.0 # how long is each window  in wpit t units
overlap_phys = 10.0 # how much they overlap

tstart_idx = int(tstart_phys / dt_phase)
tend_idx   = int(tend_phys / dt_phase)

efield_start_idx = int(tstart_phys / dt_field)
efield_end_idx = int(tend_phys / dt_field)

print(f"Phase space indices: {tstart_idx} to {tend_idx} (total: {tend_idx - tstart_idx})")
print(f"Field indices: {efield_start_idx} to {efield_end_idx} (total: {efield_end_idx - efield_start_idx})")
print(f"Ratio dt_phase/dt_field: {dt_phase/dt_field}")

window_steps = int(window_phys / dt_phase)
overlap_steps = int(overlap_phys / dt_phase)


print("tstart_idx =", tstart_idx)
print("tend_idx   =", tend_idx)
print("window_steps =", window_steps)
print("overlap_steps =", overlap_steps)
print("tend_idx - window_steps =", tend_idx - window_steps)


if tend_idx - window_steps <= tstart_idx:
    raise ValueError("No valid time windows: reduce window_phys or increase analysis interval")

time_windows = np.arange(tstart_idx,tend_idx - window_steps, overlap_steps)

#-----------------------bound frequency--------------------
time_full = np.arange(Nt_field) * DT_coeff * write_interval * mfactor
dt_snapshot = time_full[1] - time_full[0]

EF_slice = EF[int(efield_start_idx):int(efield_end_idx), :]
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
threshold = np.max(power_spectrum) * 0.01

peaks, _ = find_peaks(power_spectrum, height=threshold, distance=2)

N_modes = min(100, len(peaks))
top_order = np.argsort(power_spectrum[peaks])[-N_modes:]
top_idx = peaks[top_order]

k_dom = k_pos[top_idx] / LDi

k_avg = np.average(k_dom)
k_max = k_pos[top_idx[np.argmax(power_spectrum[top_idx])]] / LDi

E2_mean_mode = np.mean(Ek2[:, top_idx], axis=0)
k_modes_eff = k_pos[top_idx] / LDi

k_eff = np.sum(k_modes_eff * E2_mean_mode) / np.sum(E2_mean_mode)

print("k_eff =", k_eff)
print("k_dominant =", k_max)
print("k_avg:", k_avg)

# -------------------- BOUNCE FREQUENCY -----------------
Ek2_dom = Ek2[:, top_idx]
sum_Ek2_dom = np.sum(Ek2_dom, axis=1)
sum_Ek2_dom_avg = np.mean(sum_Ek2_dom)

E_eff = np.sqrt(sum_Ek2_dom_avg)
omega_b = np.sqrt((e * k_eff * E_eff) / mb)######
omega_b_norm = omega_b / wpi

print("omega_b / wpi =", omega_b_norm)
print("Tb_effective =", 2 * np.pi / omega_b_norm)

Tb_eff = 2.0 * np.pi / omega_b_norm

#--------------------------------------------------------

#------------------- STORAGE ----------------------------
bin_acf_sum   = np.zeros((nv_bins, window_steps))
bin_acf_count = np.zeros(nv_bins)

### do once 
# shape: (Nt_particle, Np, 2)
phase_all = np.stack([particle_group[k][:num_particles_to_scan] for k in phasespace_keys])

# ------------------- LOOP ---------------------------
for t0 in tqdm(time_windows, desc="Time windows"):

    t1 = t0 + window_steps

    # velocities at window start
    vel_t0 = phase_all[t0, :, 1]

    for ib in range(nv_bins):

        v_low = v_edges[ib]
        v_high = v_edges[ib + 1]

        in_bin = (vel_t0 >= v_low) & (vel_t0 < v_high)
        p_idxs = np.where(in_bin)[0]

        if p_idxs.size == 0:
            continue

        for p in p_idxs:

            xp_traj = phase_all[t0:t1, p, 0]
            Nt_w = len(xp_traj)

            E_traj = np.zeros(Nt_w)

            for j in range(Nt_w):
                f_idx = field_indices[t0 + j]
                xp = xp_traj[j]

                ix = int(np.floor(xp)) % Nx
                dxp = xp - ix
                ixn = (ix + 1) % Nx
                
                #CIC 
                E_traj[j] = (EF[f_idx, ix] * (1 - dxp) + EF[f_idx, ixn] * dxp)

            E_env = E_traj - np.mean(E_traj)

            ac = np.correlate(E_env, E_env, mode="full")[Nt_w - 1:]
            if ac[0] == 0:
                continue

            ac /= ac[0]

            bin_acf_sum[ib] += ac
            bin_acf_count[ib] += 1

# ---------------- AVERAGE ACF ----------------
acf = np.zeros_like(bin_acf_sum)
for ib in range(nv_bins):
    if bin_acf_count[ib] > 0:
        acf[ib] = bin_acf_sum[ib] / bin_acf_count[ib]

tau = np.arange(window_steps) * dt_phase

# ---------------- DECORRELATION TIME ----------------
def exp_decay(t, tc):
    return np.exp(-t / tc)

tau_c = np.zeros(nv_bins)
Ku = np.zeros(nv_bins)

for ib in range(nv_bins):
    if bin_acf_count[ib] < 10:
        continue

    fit_end = int(0.3 * window_steps)

    try:
        popt, _ = curve_fit(
            exp_decay,
            tau[:fit_end],
            acf[ib, :fit_end],
            p0=(5 * dt_phase),
            maxfev=5000
        )
        tau_c[ib] = abs(popt[0])
        Ku[ib] = tau_c[ib] / Tb_eff
    except:
        pass

# ---------------- SAVE ----------------
savefile = pjoin(path_out, f"kubo_ensemble_time_{particle_type}.h5")
with h5py.File(savefile, "w") as hf:
    hf.create_dataset("v_centers", data=v_centers)
    hf.create_dataset("tau", data=tau)
    hf.create_dataset("acf", data=acf)
    hf.create_dataset("tau_c", data=tau_c)
    hf.create_dataset("Ku", data=Ku)
    hf.create_dataset("counts", data=bin_acf_count)


# ---------------- PLOT ----------------
plt.figure(figsize=(7,5))
plt.plot(v_centers, Ku, "o-", lw=2)
plt.xlabel(r"$v$")
plt.ylabel(r"$Ku(v)$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_out, "Ku_vs_v.png"), dpi=200)
plt.show()
