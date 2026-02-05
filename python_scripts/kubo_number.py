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

#----------------constanst--------------------

eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#----------------------------------------------

if len(sys.argv) != 3:
    print("Usage: python Dth_auto_full.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]
species_name = sys.argv[2]
file_name = 'result.h5'

output_dir = 'turbulence_plots'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

#open hdf file
f = h5py.File(pjoin(path, file_name), 'r')

#---------------------Read Metadata-----------------------------------------
metadata_group = f['/metadata']

# Read Metadata
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

species_metadata_path = f'/metadata_species/{species_name}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']
species_spwt = metadata_species.attrs['spwt']

print(species_spwt)
#----------------------------------------------

#-------------------------norm------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

DT = DT_coeff * (1.0 / wpi)

L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi


efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12


# Load Electric Field Data
electric_field_data = []

time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
time_full = np.arange(Nt) * DT_coeff * write_interval*mfactor #normalized time

#normailized time between two field snapshots similar to particle pos snapshots
dt_field = DT_coeff * write_interval * mfactor

# Spatial FFT and Peak Finding
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx) #normalized k we havenot used dx = debye lenght so k is normalized
k_modes = kfft[:Nx // 2] # Positive k modes
F_k_pos = F_k[:, :Nx // 2] # Positive k modes

power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)

print("peak indices dominant k :", peaks)
print("Corresponding k values:", k_modes[peaks])
max_idx = peaks[np.argmax(power_spectrum[peaks])]

# Isolate the dominant mode's wavenumber and time-series data
k_dominant = k_modes[max_idx]
E_k_t = F_k_pos[:, max_idx]  # Amplitude E(k, t)

dt_snapshot = time_full[1] - time_full[0]
omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)

omega_pos_freqs = omega_freqs[:Nt // 2] # Positive frequencies
omega_pos_spectrum = omega_spectrum[:Nt // 2] # Positive frequencies modes
omega_peak_idx = np.argmax(np.abs(omega_pos_spectrum))
omega_r = omega_pos_freqs[omega_peak_idx]

v_phase = omega_r / k_dominant  #normalized phase velocity

print(f"Dominant Wavenumber (k): {k_dominant:.4f}")
print(f"Dominant Frequency (ω_r): {omega_r:.4f}")
print(f"Phase Velocity (v_φ = ω_r/k * vthi): {v_phase:.4f}")

# Bounce frequency calculation (normalized ω_b / ω_pe ≈ sqrt(k * |E_k|))
# E_k_amp is E_0(t) amplitude of dominant k mode

E_k_amp = np.abs(E_k_t)* efield_norm * (2/np.sqrt(NC))
print("E_k amplitude units:", E_k_amp)



"""
fft result are proportional to have to scale by 1/N or in case of ortho(where persaval idenity is preserved) 1/sqrt(N) 
and a multiplication factor of 2 comes because energy is divided into both positive and negative spectra , so positive spectra carry half and negative half.

"""

q_mass_ratio = e / species_mass
#q_mass_ratio = ion_charge / electron_mass

omega_b_avg = (np.sqrt((k_dominant / L )* np.mean(E_k_amp)*(q_mass_ratio)))/W
omega_b_max = (np.sqrt((k_dominant /L )* np.max(E_k_amp)*(q_mass_ratio)))/W
print(f"Average Bounce Frequency (ω_b / ω_pe): {omega_b_avg:.4f}")
print(f"Max Bounce Frequency (ω_b / ω_pe): {omega_b_max:.4f}")



####new
# ------------------ Lagrangian Autocorrelation ------------------
particle_start, particle_end = 5, 8 #100000, 200000

if((particle_end-particle_start) > species_num):
    print("particel number exceded total particle. error exiting ?")
    exit(-1)

# particle_group.keys() gives = ["0", "1*write_interval_phase", "2*write_interval_phase", ..]
#sort them and convert them to intergers
particle_group = f[f"particle_{species_name}"]
phasespace_data_index = sorted(particle_group.keys(), key=int)

Nt_particle = len(phasespace_data_index)

#time step between two consecutive particle snapshots, in normalized time units.(if write_int_phase is say 1000 then snapshot are taken at 0,1000,2000 
# so there exists a gap of 1000*DT_coeff*mfactor between two consecutive snapshots)
dt_particle = DT_coeff * write_interval_phase * mfactor

Nt_field = EF.shape[0]

field_indices = np.round(np.linspace(0, Nt_field - 1, Nt_particle)).astype(int)

acfs_list = []

f = h5py.File(pjoin(path, file_name), 'r')
particle_group = f[f"particle_{species_name}"]

# Avarge the ACF over multiple particles ( or if possible all the particle) to get statistically significant result

for p in tqdm(range(particle_start, particle_end + 1), desc="Computing ACF", unit="particle"):
#for p in range(particle_start, particle_end + 1):
    particle_x = np.array([particle_group[k][p, 0] for k in phasespace_data_index])
    E_traj = np.zeros(Nt_particle)

    for t in range(Nt_particle):
        idx = field_indices[t]
        xp = particle_x[t]
        i = int(np.floor(xp))
        i = np.clip(i, 0, len(x) - 2)
        dxp = xp - i
        # CIC like interpolation to calculate E at particle position
        E_traj[t] = EF[idx, i] * (1 - dxp) + EF[idx, i + 1] * dxp

    #E_env = np.abs(E_traj) - np.mean(np.abs(E_traj))
    E_env = (E_traj - np.mean(E_traj))*efield_norm
    ac = np.correlate(E_env, E_env, mode='full')[len(E_env) - 1:]
    ac /= ac[0]
    acfs_list.append(ac)

# Particle Average ACF
ac_norm = np.mean(acfs_list, axis=0)
lags_time = np.arange(Nt_particle) * dt_particle

#Integral length scale (Wiki)
if np.any(ac_norm <= 0):
    cutoff = np.where(ac_norm <= 0)[0][0] # [][] give array fo index wjich meet ocndition and and [][0] gives 1st indices
else:
    cutoff = len(ac_norm)


print("acf at cutoff", ac_norm[cutoff])
print("lags time cutt off",lags_time[cutoff])
tau_int = np.trapz(ac_norm[:cutoff], lags_time[:cutoff])
print("Integral tau:", float(tau_int))

# Exponential fit
def exp_decay(t, a, tau, c):
    return a * np.exp(-t / tau) + c

fit_end = int(0.10 * len(ac_norm))
popt, _ = curve_fit(exp_decay, lags_time[:fit_end], ac_norm[:fit_end], p0=(1, lags_time[fit_end//4], 0))
tau_exp = abs(popt[1])

#print("cutt off",lags_time[:fit_end])

# Kubo number
Ku_int = (omega_b_avg * tau_int)/(2*np.pi)
Ku_exp = (omega_b_avg * tau_exp)/(2*np.pi)

print(f"τ_int = {tau_int:.4e}, τ_exp = {tau_exp:.4e}")
print(f"Ku_int = {Ku_int:.4f}, Ku_exp = {Ku_exp:.4f}")

# Plotting
plt.figure()
plt.plot(lags_time, ac_norm, label='ACF')
# Shade the integration region
#plt.fill_between(lags_time[:cutoff],ac_norm[:cutoff],alpha=1,label="Integrated area")
#plt.axvline(tau_int, color='r', linestyle='--', label=fr"$\tau_{{int}}={tau_int:.2f}$")
plt.plot(lags_time[:fit_end], exp_decay(lags_time[:fit_end], *popt), 'g--',label=fr"Exp fit ($\tau_{{exp}}={tau_exp:.2f}$)")
plt.xlabel(r"$\omega_{pi}t$")
plt.ylabel("Autocorrelation")
plt.title(f"Particle avaerage ACF for {species_name}")
plt.legend(); 
plt.grid()
plt.tight_layout()
plt.savefig(pjoin(path_fig, f"acf_particle_avg_for10p_{species_name}.png"))
plt.show()
