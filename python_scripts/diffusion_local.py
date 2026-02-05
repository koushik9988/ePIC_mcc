"""
@kaushik kalita
@Date : 02/12/2025
Works with datasets of 1st paper/work.
"""
import numpy as np
import h5py
import sys
import os
import matplotlib.pyplot as plt
from os.path import join as pjoin
from scipy.constants import value as constants
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from mpl_toolkits.mplot3d import Axes3D


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

output_dir = 'new_analysis'
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


DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

#DT = DT_coeff * (1.0 / wpi)

L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi

DT = DT_coeff * (1.0 / W)
efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12


dt_phase = write_interval_phase * DT_coeff * mfactor

# ------------------- Load Particle Data -------------------
particle_group_name = f"particle_{particle_type}"

if particle_group_name not in f:
    print(f"Error: {particle_group_name} group not found in HDF5 file.")
    sys.exit(1)

# Get keys, convert to integers, and sort
keys = f[particle_group_name].keys()
phase_idx_list = []
for k in keys:
    phase_idx_list.append(int(k))

phase_idx = np.array(sorted(phase_idx_list))

# Load velocity data
vel_list = []
for j in phase_idx:
    # Dataset path: particle_type/timestep
    ds = f[f"{particle_group_name}/{j}"]
    # Column 1 is usually velocity
    v_data = ds[:, 1].astype(float)
    vel_list.append(v_data)

vel = np.array(vel_list)
nt, npart = vel.shape

print(f"Loaded {nt} phase-snapshots with {npart} particles each.")

# ------------------- Velocity Binning -------------------
v_flat = vel.flatten()


vmin = np.percentile(v_flat, 5)
vmax = np.percentile(v_flat, 95)
nv_bins = 50

v_edges = np.linspace(vmin, vmax, nv_bins + 1)
# Calculate centers manually for clarity
v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])


# ---------------- Time Window & FFT Calculation --------------------
electric_field_data = []


field_keys = f['fielddata/efield'].keys()
time_steps = sorted([int(k) for k in field_keys])
    
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])
    
EF = np.vstack(electric_field_data)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
dx = x[1] - x[0]
    
# Normalized time array
time_full = np.arange(Nt) * DT_coeff * write_interval * mfactor
    
# Perform FFT
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)
    
# Keep positive modes
k_modes = kfft[:Nx // 2]
F_k_pos = F_k[:, :Nx // 2]
    
# Power spectrum
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
    
# Find peaks
peak_height = np.max(power_spectrum) * 0.04
peaks, _ = find_peaks(power_spectrum, height=peak_height, distance=2)
    
if len(peaks) > 0:
    # Find the index in 'peaks' that corresponds to the max power
    peak_vals = power_spectrum[peaks]
    max_peak_loc = np.argmax(peak_vals)
    max_idx = peaks[max_peak_loc]
else:
    max_idx = np.argmax(power_spectrum)
        
k_dominant = k_modes[max_idx]
E_k_t = F_k_pos[:, max_idx]
    
    # Phase velocity calculation
dt_snapshot = time_full[1] - time_full[0]
omega_spectrum = np.fft.fft(E_k_t)
omega_freqs = 2 * np.pi * np.fft.fftfreq(Nt, dt_snapshot)
    
omega_pos_freqs = omega_freqs[:Nt // 2]
omega_pos_spectrum = omega_spectrum[:Nt // 2]
    
omega_peak_idx = np.argmax(np.abs(omega_pos_spectrum))
omega_r = omega_pos_freqs[omega_peak_idx]
    

v_phase = omega_r / k_dominant
   
    
# Determine start time based on Max Field Amplitude
Ek_amp = np.abs(E_k_t) * efield_norm * (2 / np.sqrt(Nx))
Ek_ampmax_idx = np.argmax(Ek_amp)
tstart_idx = time_full[Ek_ampmax_idx]

# ---------------- Define Analysis Time Window ----------------
field_slice_for_lenght_diffusion = 30

t0_idx = int(round(tstart_idx / dt_phase))
n_ti = int(round((tstart_idx + field_slice_for_lenght_diffusion) / dt_phase))

t0_idx = int(round(100/ dt_phase))
n_ti = int(round((490)/ dt_phase))

print(f"Diffusion Analysis Time Window: Start Index = {t0_idx}")
print(f"Diffusion Analysis Time Window: End Index = {n_ti}")

# Create range of time indices
ti_idxs_raw = np.arange(t0_idx, t0_idx + n_ti)
# Filter indices to ensure they don't exceed data length
ti_idxs = []
for t in ti_idxs_raw:
    if t < nt:
        ti_idxs.append(t)
ti_idxs = np.array(ti_idxs)

# Define Lag Steps (Tau)
part1 = np.arange(1, 6)
part2 = np.arange(6, 31, 2)
part3 = np.arange(32, 101, 8)
tau_steps = np.unique(np.concatenate([part1, part2, part3]))
tau_phys = tau_steps * dt_phase

# ---------------- Initialize Storage Arrays ----------------
D1 = np.zeros(nv_bins)
D2 = np.zeros(nv_bins)
D3 = np.zeros(nv_bins)
D2_from_sigma = np.zeros(nv_bins)
counts = np.zeros(nv_bins, int)

n_pdf_bins = 100
all_var = []

#----
#all_pdfs[ib, i_tau, :]  → PDF curve of Δv distribution for bin ib at lag tau_steps[i_tau]
all_pdfs = np.zeros((nv_bins, len(tau_steps), n_pdf_bins))
#all_sigma_fit[ib, i_tau] = σ(τ)
all_sigma_fit = np.zeros((nv_bins, len(tau_steps)))
all_mu_fit = np.zeros((nv_bins, len(tau_steps)))
all_amp_fit = np.zeros((nv_bins, len(tau_steps))) 
all_D2_from_sigma_tau = np.zeros((nv_bins, len(tau_steps)))
all_sigma_msd = np.zeros((nv_bins, len(tau_steps)))  ##store sigma from MSD for all v-bins


# Store the exact x-axis (centers) for every single bin to prevent shift errors
pdf_centers_storage = [] 

# Number of points to use for linear slope fitting
Nfit = len(tau_steps) // 4
if Nfit < 6:
    Nfit = 6

#Gaussian function for curve_fit.
def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

# ================= MAIN LOOP OVER VELOCITY BINS =================
for ib in range(nv_bins):
    v_low = v_edges[ib]
    v_high = v_edges[ib + 1]

    # Dictionary to hold velocity differences for each tau
    dv_per_tau = {}
    for tau in tau_steps:
        dv_per_tau[tau] = []


    for ti in ti_idxs:
        if ti >= nt:
            continue
            
        # Identify particles in this bin at time ti
        in_bin_mask = (vel[ti] >= v_low) & (vel[ti] < v_high)
        idxs = np.flatnonzero(in_bin_mask)
        
        if idxs.size == 0:
            continue

        for tau in tau_steps:
            t2 = ti + tau
            if t2 >= nt:
                continue
            
            # Calculate dv = v(t + tau) - v(t)
            v_final = vel[t2, idxs]
            v_initial = vel[ti, idxs]
            dv = v_final - v_initial
            
            dv_per_tau[tau].append(dv)

    #Flatten data for range determination
    all_dv_samples_for_bin = []
    for tau in tau_steps:
        if len(dv_per_tau[tau]) > 0:
            concatenated_dv = np.concatenate(dv_per_tau[tau])
            all_dv_samples_for_bin.append(concatenated_dv)

    all_dv_concat = np.concatenate(all_dv_samples_for_bin)
    
    #Establish PDF Bins for THIS SPECIFIC velocity slice
    dv_min = np.min(all_dv_concat)
    dv_max = np.max(all_dv_concat)
    dv_range = dv_max - dv_min
    
   
    buffer_width = 0.05 * dv_range
    
    pdf_bins = np.linspace(dv_min - buffer_width, dv_max + buffer_width, n_pdf_bins + 1)
    
    # Calculate centers
    pdf_bin_centers = 0.5 * (pdf_bins[:-1] + pdf_bins[1:])
    
    # IMPORTANT: Store these centers to align plots later
    pdf_centers_storage.append(pdf_bin_centers)

    msd_list = []
    var_list = []
    
    # Process each Tau step
    for i_tau, tau in enumerate(tau_steps):
        dv_list = dv_per_tau[tau]
                 
        dv_cat = np.concatenate(dv_list)
        
        # Calculate MSD and VAR
        msd_val = np.mean(dv_cat**2)
        var_val = np.var(dv_cat)
        
        msd_list.append(msd_val)
        var_list.append(var_val)
        
        # Calculate Histogram
        hist_vals, _ = np.histogram(dv_cat, bins=pdf_bins, density=True)
        all_pdfs[ib, i_tau, :] = hist_vals
        
        # Find the Peak (Mode) index
        peak_idx = np.argmax(hist_vals)
        
        # This fixes the shifting issue caused by heavy tails.
        mu0 = pdf_bin_centers[peak_idx] 
        A0 = hist_vals[peak_idx]
        sigma0 = np.std(dv_cat) 
        
        p0 = [A0, mu0, sigma0]

        #Define Bounds
        lower_bounds = [0, pdf_bins[0], 1e-12]
        upper_bounds = [np.inf, pdf_bins[-1], np.inf]
        bounds_tuple = (lower_bounds, upper_bounds)

        

#*********************************************************************************************************************

        #Perform Fit
        #Curve_fit(f, xdata, ydata, p0=None, sigma=None, absolute_sigma=False, check_finite=None, 
        # bounds=(-inf, inf), method=None, jac=None, *, full_output=False, nan_policy=None, **kwargs)
        ###(yse try block such if one Gaussian fit fails, it does not stop the entire process)
        try:
            popt, pcov = curve_fit(gauss, pdf_bin_centers, hist_vals,p0=p0, bounds=bounds_tuple, maxfev=8000)

            A_fit, mu_f, sigma_f = popt
            all_amp_fit[ib, i_tau]    = A_fit
            all_mu_fit[ib, i_tau]     = mu_f
            all_sigma_fit[ib, i_tau]  = abs(sigma_f)
            tau_val = tau_phys[i_tau]
            all_D2_from_sigma_tau[ib, i_tau] = (sigma_f**2) / (2.0 * tau_val)

        except RuntimeError:
        # Gaussian fit failed → non-Gaussian distribution
            all_amp_fit[ib, i_tau]    = np.nan
            all_mu_fit[ib, i_tau]     = np.nan
            all_sigma_fit[ib, i_tau]  = np.nan
            all_D2_from_sigma_tau[ib, i_tau] = np.nan
            continue   # move to next tau-step without giving errro
        ###

        """
        popt, pcov = curve_fit(gauss, pdf_bin_centers, hist_vals, p0=p0, bounds=bounds_tuple, maxfev=5000)   
        A_fit, mu_f, sigma_f = popt
        # Store results
        all_amp_fit[ib, i_tau] = A_fit
        all_mu_fit[ib, i_tau] = mu_f
        all_sigma_fit[ib, i_tau] = abs(sigma_f)
        
        # Calculate D2 based on sigma
        tau_val = tau_phys[i_tau]
        all_D2_from_sigma_tau[ib, i_tau] = (sigma_f**2) / (2.0 * tau_val)
        """
#*********************************************************************************************************************

    #Store Variance Array
    all_var.append(np.array(var_list))

    msd_array = np.array(msd_list)
    var_array = np.array(var_list)

    sigma_msd = np.sqrt(msd_array)

    all_sigma_msd[ib, :] = sigma_msd  ##save sigma(MSD) vs tau for this velocity bin

    # 1) Diffusion from VAR slope:  Var = 2Dτ + C
    slope_var, intercept_var = np.polyfit(tau_phys[:Nfit], var_array[:Nfit], 1)
    D_var = slope_var / 2.0
    D1[ib] = D_var 


    # 2) Diffusion from MSD:  MSD = 2Dτ + C
    slope_msd, intercept_msd = np.polyfit(tau_phys[:Nfit], msd_array[:Nfit], 1)
    D_msd = slope_msd / 2.0
    D2[ib] = D_msd 


    # 3) Diffusion from sigma² using curve-fit: σ² = 2Dτ + C
    sigma_fit_arr = all_sigma_fit[ib, :Nfit]
    sigma2 = sigma_fit_arr**2
    slope_sigma, intercept_sigma = np.polyfit(tau_phys[:Nfit], sigma2[:Nfit], 1)
    D_sigma = slope_sigma / 2.0
    D3[ib] = D_sigma
    D2_from_sigma[ib] = D_sigma

    if ib == 2:

        plt.figure()

        # Scatter values
        plt.scatter(tau_phys[:Nfit], msd_array[:Nfit], s=35, color='blue', label = r'$\sigma_{MSD}$')
        plt.scatter(tau_phys[:Nfit], sigma2[:Nfit],s=35, color='red', label = r'$\sigma_{Gauss}$')

        # ======== Fitted straight lines ========
        sigma2_msd_fit   = slope_msd*tau_phys[:Nfit] + intercept_msd
        sigma2_gauss_fit = slope_sigma*tau_phys[:Nfit] + intercept_sigma

        plt.plot(tau_phys[:Nfit], sigma2_msd_fit, color='blue', linewidth=2)

        plt.plot(tau_phys[:Nfit], sigma2_gauss_fit, color='red', linewidth=2)


        plt.xlabel(r'$\tau$', fontsize=14)
        plt.ylabel(r'$\sigma(\Delta v)$', fontsize=14)
       
        plt.grid(alpha=0.35)
        plt.legend(fontsize=11)
        plt.tight_layout()
        plt.show()
        plt.savefig(pjoin(path_fig, f"sigma_vs_tau{ib}.png"), dpi=900)


    # Count total samples used
    #total_samples = 0
    #for t in tau_steps:
    #    total_samples += len(dv_per_tau[t])
    #counts[ib] = total_samples



# ----------------- PLOTTING -----------------
#Diffusion Coefficients vs Velocity
plt.figure()
plt.plot(v_centers, D1, '-o', label=r'$D_1 = \langle (\Delta v)^2 \rangle / (2\tau)$')
#plt.plot(v_centers, D2, '-s', label=r'$D_2 = \mathrm{Var}(\Delta v)/(2\tau)$')
plt.plot(v_centers, D2_from_sigma, '-^', label=r'$D_2 = \sigma / 2\tau$')
plt.xlabel(r'$v$')
plt.ylabel(r'$D(v)$')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
save_path_1 = pjoin(path_fig, f"D_local_vs_v_{particle_type}.png")
plt.savefig(save_path_1, dpi=180)
plt.show()

#Variance Evolution for all bins
#plt.figure()
#for ib in range(len(all_var)):
#    current_var = all_var[ib]
#    # Only plot if we have a full dataset
#    if len(current_var) == len(tau_phys):
#        plt.plot(tau_phys, current_var, alpha=0.35, linewidth=1)

#plt.xlabel(r"$\tau$")
#plt.ylabel(r"$\delta v^2$ (Var)")
#plt.title("VAR(τ) vs τ")
#plt.grid(True, alpha=0.3)
#lt.tight_layout()
#save_path_2 = pjoin(path_fig, "VAR_vs_tau_ALL_bins.png")
#lt.savefig(save_path_2, dpi=160)
#plt.show()

#f.close()

# Select the velocity bin to plot
# Priority: Phase velocity -> Center of bins

#extra plots for PDF and sigma vs tau

if np.abs(v_phase) > 0:
    v_choice = v_phase
else:
    mid_index = nv_bins // 2
    v_choice = v_centers[mid_index]

v_choice = v_phase
dist = np.abs(v_centers - v_choice)
sample_ib = np.argmin(dist)
v0_center = v_centers[sample_ib]

print(f"Plotting PDF for bin index: {sample_ib}, v0={v0_center:.3f}")


dv_centers = pdf_centers_storage[sample_ib] 
tau_indices = [1, len(tau_steps)//2, max(1, len(tau_steps)-2)]
tau_indices = [1,2,3,4,5,6,7,8,9,10]#,11,12,13,14,15,16,17,18,19,20]
#tau_indices = [1,5,10,20]

# Plot PDF + Gaussian Fit
plt.figure()

for idx in tau_indices:
    # Get Histogram Data
    pdf_vals = all_pdfs[sample_ib, idx, :]
    
    # Get Fit Parameters
    sigma_f = all_sigma_fit[sample_ib, idx]
    mu_f = all_mu_fit[sample_ib, idx]
    A_f = all_amp_fit[sample_ib, idx]

    tau_val = tau_steps[idx] * dt_phase

    # Plot Data Dots
    plt.plot(dv_centers, pdf_vals, '-', markersize=4, alpha=0.6,label=rf' $\tau$ = {idx}')

    x_smooth_min = dv_centers[0]
    x_smooth_max = dv_centers[-1]
    x_smooth = np.linspace(x_smooth_min, x_smooth_max, 200)
        
    fit_y = gauss(x_smooth, A_f, mu_f, sigma_f)
        
    #plt.plot(x_smooth, fit_y, 'o', linewidth=1.5,label=rf'Fit $\mu$={mu_f:.2f}')

plt.xlabel(r'$\Delta v$')
plt.ylabel(r'$P_{v_0}(\Delta v, \tau)$')
#plt.yscale('log')
plt.title(rf'PDF Fits  for $v_0\approx{v0_center:.3f}$')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
save_path_3 = pjoin(path_fig, f"PDF_with_Gaussian_fit_v0_bin_{sample_ib}.png")
plt.savefig(save_path_3, dpi=180)
plt.show()

# Plot Sigma vs Tau
plt.figure()
sigma_tau = all_sigma_fit[sample_ib, :]

# Filter zero/invalid sigmas
#mask_sigma = sigma_tau > 1e-12
plt.plot(tau_phys, sigma_tau, 'o-', linewidth=2)

plt.xlabel(r'$\tau$')
plt.ylabel(r'$\sigma(\Delta v)$')
plt.title(rf'sigma(τ) evolution')
plt.grid(True, alpha=0.3)
plt.tight_layout()
save_path_4 = pjoin(path_fig, f"sigma_from_P_vs_tau_bin_{sample_ib}.png")
plt.savefig(save_path_4, dpi=180)
plt.show()



pdf_plot = all_pdfs[sample_ib,:,:]

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection='3d')

num_skip = 5  # plot every 5th tau

pdf_plot = pdf_plot / np.max(pdf_plot)

for i in range(0,len(tau_phys), num_skip):
    #ax.plot(dv_centers, pdf_plot[i,:], tau_phys[i] * np.ones_like(dv_centers),
            #lw=1.2)
    ax.plot(dv_centers, tau_phys[i] * np.ones_like(dv_centers), pdf_plot[i,:],
            lw=1.8)
ax.grid(False)
ax.set_xlabel(r'$\Delta v$')
ax.set_zlabel(r'$ P(\Delta v,\tau)$')
ax.set_ylabel(r'$\tau$')
plt.tight_layout()
plt.show()
plt.savefig(pjoin(path_fig, f"PDF_v0_bin_{sample_ib}.png"), dpi=800)



###save all
save_h5 = pjoin(path_fig, f"local_velocity_diffusion_{particle_type}.h5")

with h5py.File(save_h5, "w") as hf:

    # ---------------- Metadata ----------------
    gmeta = hf.create_group("metadata")
    gmeta.attrs["particle_type"] = particle_type
    gmeta.attrs["nv_bins"] = nv_bins
    gmeta.attrs["n_tau"] = len(tau_steps)
    gmeta.attrs["dt_phase"] = dt_phase
    gmeta.attrs["v_phase"] = v_phase
    gmeta.attrs["analysis_window_start"] = t0_idx
    gmeta.attrs["analysis_window_end"] = n_ti
    gmeta.attrs["Nfit"] = Nfit
    gmeta.attrs["DT_phase"] = dt_phase

    # ---------------- Axes info for ploting later ----------------
    axes = hf.create_group("axes")
    axes.create_dataset("v_edges", data=v_edges)
    axes.create_dataset("v_centers", data=v_centers)
    axes.create_dataset("tau_steps", data=tau_steps)
    axes.create_dataset("tau_phys", data=tau_phys)

    # Each velocity bin has its own Δv axis
    pdf_axes = axes.create_group("pdf_centers")
    for ib in range(nv_bins):
        pdf_axes.create_dataset(f"bin_{ib}", data=pdf_centers_storage[ib])

    # ---------------- PDFs ----------------
    gpdf = hf.create_group("pdfs")
    gpdf.create_dataset("P_dv_tau_v", data=all_pdfs, compression="gzip")

    # ---------------- Gaussian Fits ----------------
    gfit = hf.create_group("gaussian_fit")
    gfit.create_dataset("A", data=all_amp_fit)
    gfit.create_dataset("mu", data=all_mu_fit)
    gfit.create_dataset("sigma", data=all_sigma_fit)

    # ---------------- Diffusion from sigma(τ) ----------------
    gDsig = hf.create_group("diffusion_from_sigma_tau")
    gDsig.create_dataset("D_tau", data=all_D2_from_sigma_tau)

    # ---------------- Final Local Diffusion Coefficients ----------------
    gD = hf.create_group("local_diffusion")
    gD.create_dataset("D_from_var", data=D1)
    gD.create_dataset("D_from_msd", data=D2)
    gD.create_dataset("D_from_sigma", data=D2_from_sigma)

    # ---------------- Variance / MSD ----------------
    gstats = hf.create_group("statistics")
    gstats.create_dataset("variance", data=np.array(all_var))
    gstats.create_dataset("sigma_from_msd", data=all_sigma_msd) 



print(f"All analysis saved to:\n{save_h5}\n")
