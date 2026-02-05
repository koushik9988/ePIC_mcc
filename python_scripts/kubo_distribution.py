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
particle_type = sys.argv[2]
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

species_metadata_path = f'/metadata_species/{particle_type}'
metadata_species = f[species_metadata_path]

species_mass = metadata_species.attrs['mass']
species_charge = metadata_species.attrs['charge']
species_temp = metadata_species.attrs['temperature']
species_spwt = metadata_species.attrs['spwt']
species_num = metadata_species.attrs['num_particles']
species_spwt = metadata_species.attrs['spwt']


q_mass_ratio = e / species_mass

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



#-----Load Electric Field Data and stackto form a 2D matrix--------------------------------------
electric_field_data = []

time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]


#--------------time array normalized time--------------------------------
time_full = np.arange(Nt) * DT_coeff * write_interval*mfactor 

#---------normailized time between two field snapshots similar to particle pos snapshots----------
dt_field = DT_coeff * write_interval * mfactor


#---------------------------------FFT-----------------------------------------------------
#----Spatial FFT and Peak Finding
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx) # Normalized k we havenot used dx = debye lenght so k is normalized
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

v_phase = omega_r / k_dominant  # Normalized phase velocity

print(f"Dominant Wavenumber (k): {k_dominant:.4f}")
print(f"Dominant Frequency (ω_r): {omega_r:.4f}")
print(f"Phase Velocity (v_φ = ω_r/k * vthi): {v_phase:.4f}")

# Bounce frequency calculation (normalized ω_b / ω_pe ≈ sqrt(k * |E_k|))
# E_k_amp is E_0(t) amplitude of dominant k mode
"""
fft result are proportional to have to scale by 1/N or in case of ortho(where persaval idenity is preserved) 1/sqrt(N) 
and a multiplication factor of 2 comes because energy is divided into both positive and negative spectra , so positive spectra carry half and negative half.

"""
E_k_amp = np.abs(E_k_t)* efield_norm * (2/np.sqrt(Nx))
print("E_k amplitude units:", E_k_amp)

#---------------trapped-time-calculation -------------------------
dEk_dt = np.gradient(E_k_t, dt_field, axis=0)
dEkmax_idx = np.argmax(dEk_dt)
buffer = 30
tstart_idx = time_full[dEkmax_idx] + buffer

print("field frame index",int(round(tstart_idx / dt_field)))
start_time = int(round(tstart_idx / dt_field)) #120 # starting time for stat of tubulent rehion in wpit units
end_time = int(round((Nt - 2) / dt_field))

start_time = max(0, start_time)
end_time = min(Nt, end_time)
if start_time >= end_time:
    raise ValueError("Invalid time window for averaging: start_time >= end_time")

sum_num = 0.0
count = 0
print("field snapshot start :",start_time,"and end time :",end_time)
for time in range(start_time, end_time):
    Ek = np.fft.fft(EF[time:time+1, :], axis=1, norm='ortho') *(2/np.sqrt(Nx))  #shape (1, 1024)
    k = 2.0 * np.pi * np.fft.fftfreq(Nx, dx)
    Ek = Ek[:, 1:Nx // 2]       
    k_vals = k[1:Nx // 2]/LDi
    Ek2 = np.abs(Ek)**2          
    sum_num += np.sum((k_vals**2) * Ek2)
    count += 1
    
avg_k2E2 = (sum_num)/count 
kE_eff = np.sqrt(avg_k2E2)

kE_eff_phys = kE_eff * efield_norm
omega_b_eff = np.sqrt((q_mass_ratio) * kE_eff_phys) / wpi
Tb_eff = 2.0 * np.pi / omega_b_eff

print(f"Bounce frequency ω_b_eff / wpi = {omega_b_eff}")

#---------------------   end-----------------------------------------

# ------------------ Lagrangian Autocorrelation ------------------
particle_start, particle_end = 0,100#499999  #inclusive

if (particle_end - particle_start + 1) > species_num:
    print("particle number exceeded total particle. error!!")
    exit(-1)

particle_group = f[f"particle_{particle_type}"]
phasespace_data_index = sorted(particle_group.keys(), key=int)
Nt_particle = len(phasespace_data_index)
dt_particle = DT_coeff * write_interval_phase * mfactor
Nt_field = EF.shape[0]
field_indices = np.round(np.linspace(0, Nt_field - 1, Nt_particle)).astype(int)

acfs_list = []

# ------------------ time window linear and non linear regime classification-----------------
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
E_k_amp_time = E_k_amp       
field_idx_peak = np.argmax(E_k_amp_time) #index in EF/time_full where E_k is maximum

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = 0.1 * E_k_amp_time[field_idx_peak] # drop below 30 percent of peak
valid_indices = np.where(E_k_amp_time[field_idx_peak:] < threshold)[0]

if len(valid_indices) > 0:
    print("Found a point where ampltude drops 10 percent of peak")
    field_idx_end_rel = valid_indices[0] 
    field_idx_end = field_idx_peak + field_idx_end_rel
else:
    # If never drops that low, take the whole series
    field_idx_end = len(E_k_amp_time) - 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_peak_field = time_full[field_idx_peak]
t_end_field = time_full[field_idx_end] 

print(f"Dominant mode peaks at field index {field_idx_peak}, time {t_peak_field:.6g}")

# Choose linear window = [0, t_peak_field]
t1 = 0.0
t2 = tstart_idx #120#t_peak_field
t3 = t_end_field

#buffer
t1_buffer = 0.0 * dt_particle          
t2_buffer = 0.0 * dt_particle
t3_buffer = 0.0* dt_particle         

t1 = max(0.0, t1 + t1_buffer)
t2 = min(time_full[-1], t2 + t2_buffer)
t3 = min(time_full[-1],t3 + t3_buffer)

######################

"""
time_full = (0,1*write_interval,2*write_interval,3*write_interval.....)*mfactor*DT_coeff =>field data snapshot
time_phase = (0,1*write_interval_phase,2*write_interval_phase.....)*mfactor*DT_coeff

extract phase space snapshot from fild snapshot => roundoff_to_integer(timefull/(write_interval_phase*mfactor*DT__coeff))

eg time_field = (0,1,2,3,4,5,6...)*100(write_int_field)*(DT_coeff*mafactor) = (0,100,200,300,400,500,600.......)
   time_phase_snapshot = (0,1,2,3,4,5,6...)*200(write_int_phase) = (0,200,400,600,800,1000,1200.....)

   phase_snapshot = roundoff(time_field/(200*DT_coeff*mfactor))

"""

#linear window
#i1 = int(round(t1 / dt_particle))
#i2 = int(round(t2 / dt_particle))

# Nonlinear window (from peak to end)
# nonlinear start = end of linear region
i1 = int(round(t2 / dt_particle))
i2 = Nt_particle-1#int(round(t3 / dt_particle)) #Nt_particle-1   #max is 300    # nonlinear end = last phase snapshot

print(f"Auto selected linear/nonlinear window: t1={t1:.6g}, t2={t2:.6g}")
print(f"Corresponding particle snapshot indices: i1={i1}, i2={i2}")

# Sanity check
if i1 < 0:
    i1 = 0
if i2 > Nt_particle:
    i2 = Nt_particle
if i1 >= i2:
    raise RuntimeError("ERRROR !!!!!!!")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Precompute time lags for this window
lags_time_window = np.arange(i2 - i1) * dt_particle

tau_int_list = []
tau_exp_list = []
Ku_list = []
velocity_mean_list = []

for p in tqdm(range(particle_start, particle_end), desc="Computing ACF", unit="particle"):
    particle_x_full = np.array([particle_group[k][p, 0] for k in phasespace_data_index])
    particle_x = particle_x_full[i1:i2]
    Nt_w = len(particle_x)


    ###new
    # Load Velocity Data 
    particle_v_full = np.array([particle_group[k][p, 1] for k in phasespace_data_index])
    particle_v_window = particle_v_full[i1:i2] 
    
    # Calculate Mean Velocity and store it
    v_mean = np.mean(particle_v_window)
    velocity_mean_list.append(v_mean) 
    ##new

    # Interpolate E-field on particle trajectory
    E_traj = np.zeros(Nt_w)
    for j in range(Nt_w):
        t_idx_particle = i1 + j
        field_idx = field_indices[t_idx_particle]

        xp = particle_x[j]
        i = int(np.floor(xp))
        #i = np.clip(i, 0, Nx - 2)no need bound
        dxp = xp - i

        E_traj[j] = EF[field_idx, i] * (1 - dxp) + EF[field_idx, i + 1] * dxp

    # Remove mean 
    E_env = (E_traj - np.mean(E_traj)) * efield_norm
    #E_env = (E_traj) * efield_norm

    # ACF in the window only
    ac = np.correlate(E_env, E_env, mode='full')[Nt_w-1:]  #or [len(E_env) - 1:]
    if ac[0] == 0: ac[0] = 1e-20 # Avoid encounter by zero error
    ac = ac / ac[0]

    # Determine where ACF crosses 0
    #Integral length scale (Wiki)
    if np.any(ac <= 0):
        cutoff_i = np.where(ac <= 0)[0][0] # [][] give array fo index wjich meet ocndition and and [][0] gives 1st indices
    else:
        cutoff_i = len(ac)

    tau_int = np.trapz(ac[:cutoff_i],lags_time_window[:cutoff_i])

    tau_int_list.append(tau_int)

    # Exponential fit
    fit_end  = int(0.10 * len(ac))
    #fit_end = max(3, int(0.15 * Nt_w))

    def exp_decay(t, a, tau, c):
        return a * np.exp(-t / tau) + c

    popt, _ = curve_fit(exp_decay,lags_time_window[:fit_end],ac[:fit_end],p0=(1.0, dt_particle, 0.0),maxfev=5000) 

    tau_exp = abs(popt[1])
    tau_exp_list.append(tau_exp)

    # Kubo number using bounce frequency
    #Ku_list.append(tau_exp * omega_b_max / (2 * np.pi))
    Ku_list.append(tau_exp / Tb_eff)

# Convert lists to arrays
tau_int_list = np.array(tau_int_list)
tau_exp_list = np.array(tau_exp_list)
Ku_list = np.array(Ku_list)
velocity_mean_list = np.array(velocity_mean_list)


# ---------------- Save per-particle results to a TXT file ----------------
txt_out = pjoin(path_fig, f"per_particle_kubo_{particle_type}.txt")

ftxt = open(txt_out, 'w')

ftxt.write("# particle_index  tau_int tau_exp Kubo <v>  Tb_eff \n")

data_iterator = zip(tau_int_list, tau_exp_list, Ku_list, velocity_mean_list)

for idx, (ti, te, ku, vm) in enumerate(data_iterator, start=particle_start):
    #line = (f"{idx}" f"{ti:.8e}" f"{te:.8e}" f"{ku:.8e}" f"{vm:.8e}"  f"{omega_b_max:.8e}"  f"{omega_b_avg:.8e}\n")
    line = f"{idx}, {ti:.8e}, {te:.8e}, {ku:.8e}, {vm:.8e}, {Tb_eff:.8e}\n"

    ftxt.write(line)

ftxt.close()
print("Saved per-particle Kubo TXT at:", txt_out)


# plot Kubo distribution 
Ku_plot = Ku_list
plt.figure(1)
plt.hist(Ku_plot, bins=100, edgecolor='black')
plt.xlabel("$Ku$")
plt.ylabel("Count")
plt.title(f"Kubo Number Distribution ({particle_type})")
kubo_png = pjoin(path_fig, f"kubo_distribution_{particle_type}.png")
plt.savefig(kubo_png, dpi=200)
print("Saved Kubo distribution plot at:", kubo_png)
plt.show()




""""


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

E_k_amp = np.abs(E_k_t)* efield_norm * (2/np.sqrt(Nx))
print("E_k amplitude units:", E_k_amp)


fft result are proportional to have to scale by 1/N or in case of ortho(where persaval idenity is preserved) 1/sqrt(N) 
and a multiplication factor of 2 comes because energy is divided into both positive and negative spectra , so positive spectra carry half and negative half.



q_mass_ratio = e / species_mass
#q_mass_ratio = ion_charge / electron_mass

omega_b_avg = (np.sqrt((k_dominant / L )* np.mean(E_k_amp)*(q_mass_ratio)))/W
omega_b_max = (np.sqrt((k_dominant /L )* np.max(E_k_amp)*(q_mass_ratio)))/W
print(f"Average Bounce Frequency (ω_b): {omega_b_avg:.4f}")
print(f"Max Bounce Frequency (ω_b): {omega_b_max:.4f}")



# ------------------ Lagrangian Autocorrelation ------------------
particle_start, particle_end = 0,19000  # inclusive

if (particle_end - particle_start) > species_num:
    print("particle number exceeded total particle. error!!")
    exit(-1)

particle_group = f[f"particle_{species_name}"]
phasespace_data_index = sorted(particle_group.keys(), key=int)
Nt_particle = len(phasespace_data_index)
dt_particle = DT_coeff * write_interval_phase * mfactor
Nt_field = EF.shape[0]
field_indices = np.round(np.linspace(0, Nt_field - 1, Nt_particle)).astype(int)

acfs_list = []

#for p in range(particle_start, particle_end + 1):
for p in tqdm(range(particle_start, particle_end), desc="Computing ACF", unit="particle"):
    particle_x = np.array([particle_group[k][p, 0] for k in phasespace_data_index])
    E_traj = np.zeros(Nt_particle)
    for t in range(Nt_particle):
        idx = field_indices[t]
        xp = particle_x[t]
        i = int(np.floor(xp))
        i = np.clip(i, 0, len(x) - 2)
        dxp = xp - i
        E_traj[t] = EF[idx, i] * (1 - dxp) + EF[idx, i + 1] * dxp

    E_env = (E_traj - np.mean(E_traj)) * efield_norm
    ac = np.correlate(E_env, E_env, mode='full')[len(E_env)-1:]
    if ac[0] == 0:
        ac[0] = 1e-20
    ac = ac / ac[0]
    acfs_list.append(ac)

acfs_arr = np.array(acfs_list)           # shape (n_particles, Nt_particle)
lags_time = np.arange(Nt_particle) * dt_particle

#helper
def safe_cutoff(ac_arr):
    if np.any(ac_arr <= 0):
        return np.where(ac_arr <= 0)[0][0]
    return len(ac_arr)

# exponential decay function
def exp_decay(t, a, tau, c):
    return a * np.exp(-t / tau) + c

# compute per-particle stats
tau_int_list = []
tau_exp_list = []
Ku_list = []

for ac in acfs_arr:
    cutoff_i = safe_cutoff(ac)
    tau_int_i = 0.0 if cutoff_i <= 1 else np.trapz(ac[:cutoff_i], lags_time[:cutoff_i])
    tau_int_list.append(tau_int_i)

    fit_end_i = max(3, int(0.10 * len(ac)))
    popt_i, _ = curve_fit(exp_decay,lags_time[:fit_end_i],ac[:fit_end_i],p0=(1.0, max(1e-8, lags_time[fit_end_i//3]), 0.0),maxfev=5000)
    tau_exp_i = abs(popt_i[1])
    tau_exp_list.append(tau_exp_i)

    Ku_i = tau_int_i * omega_b_avg / (2 * np.pi)
    ku_exp = tau_exp_i * omega_b_avg / (2 * np.pi)

    #Ku_list.append(Ku_i)
    Ku_list.append(ku_exp)

tau_int_list = np.array(tau_int_list)
tau_exp_list = np.array(tau_exp_list)
Ku_list = np.array(Ku_list)

# save CSV
import csv
csv_out = pjoin(path_fig, f"per_particle_kubo_{species_name}.csv")
with open(csv_out, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['particle_index', 'tau_int', 'tau_exp', 'kubo'])
    for idx, (ti, te, ku) in enumerate(zip(tau_int_list, tau_exp_list, Ku_list), start=particle_start):
        writer.writerow([idx, float(ti), float(te), float(ku)])
print("Saved per-particle Kubo CSV at:", csv_out)

# plot Kubo distribution 
Ku_plot = Ku_list
plt.figure(1)
plt.hist(Ku_plot, bins=100, edgecolor='black')
plt.xlabel("$Ku$")
plt.ylabel("Count")
plt.title(f"Kubo Number Distribution ({species_name})")
kubo_png = pjoin(path_fig, f"kubo_distribution_{species_name}.png")
plt.savefig(kubo_png, dpi=200)
print("Saved Kubo distribution plot at:", kubo_png)
plt.show()


#plt.figure(2)
#plt.hist(Ku_plot, bins=100, edgecolor='black')
#plt.semilogy()
#plt.xlabel("$Ku$")
#plt.ylabel("Count")
#kubo_png = pjoin(path_fig, f"kubo_distribution_log_{species_name}.png")
#plt.savefig(kubo_png, dpi=200)
#print("Saved Kubo distribution plot at:", kubo_png)
#plt.show()
"""