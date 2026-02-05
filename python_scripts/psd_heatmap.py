import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os

from scipy.interpolate import interp1d


if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data_directory>")
    sys.exit(1)

path = sys.argv[1]
file_name = 'result.h5'

plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

f = h5py.File(pjoin(path, file_name), 'r')

metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']

iontemp = metadata_ion.attrs['temperature']
spwt_ion = metadata_ion.attrs['spwt']
ion_num = metadata_ion.attrs['num_particles']
ion_mass = metadata_ion.attrs['mass']
ion_charge = metadata_ion.attrs['charge']

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


L = LDe

if normscheme == 2:
    L = LDi

mfactor = wpi/wpe 
if normscheme in [1, 2]:
    mfactor = 1

print("mfactor:", mfactor)
print("normalization scheme:", normscheme)



# E-field normalization
efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12

#load field
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for ts in time_steps:
    EF_data = f[f'fielddata/efield/{ts}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx, endpoint=False)
dx = x[1] - x[0]

print(f"Data shape: EF({Nt}, {Nx})")

# Spatial FFT
F_k = np.fft.fft(EF, axis=1, norm='ortho') * (efield_norm)*(2/np.sqrt(NC)) # sclae amplitude
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)
k_pos = kfft[:Nx // 2]
F_k_pos = F_k[:, :Nx // 2]


# Temporal FFT
F_k_t = np.fft.fft(F_k_pos, axis=0, norm='ortho')
omega_fft = 2 * np.pi * np.fft.fftfreq(Nt, DT_coeff * write_interval * mfactor)
omega_pos = omega_fft[:Nt // 2]
F_k_t_pos = F_k_t[:Nt // 2, :] #amplitide

# SPECTRUM (ω–k) 
dispersion_power = np.abs(F_k_t_pos)**2
dispersion_power /= np.max(dispersion_power)  # normalize

####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Find dominant (w-k) pair
# FIND DOMINANT k 
psd = np.mean(np.abs(F_k_pos)**2, axis=0) # axis = 0 mean over time
peaks_k, _ = find_peaks(psd, height=np.max(psd) * 0.01, distance=2)
dominant_k = k_pos[peaks_k]

#FIND DOMINANT ω FOR EACH k 
dominant_omega = []

for idx in peaks_k:
    #Get PSD vs time for each k
    psd_vs_omega = np.abs(F_k_t_pos[:, idx])**2
    """
    The colon : means take all rows 
    idx selects a specific column → a particular k value.
    """
    peak_omega_idx, _ = find_peaks(psd_vs_omega, height=np.max(psd_vs_omega) * 0.01, distance=1)
    
    # Take the strongest peak in omega for this k
    if len(peak_omega_idx) > 0:
        strongest_idx = peak_omega_idx[np.argmax(psd_vs_omega[peak_omega_idx])]
        dominant_omega.append(omega_pos[strongest_idx])
    else:
        dominant_omega.append(np.nan)  # in case no peak is found

dominant_omega = np.array(dominant_omega)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(dominant_k)

fig, ax = plt.subplots()
extent = [k_pos[0], k_pos[-1], omega_pos[0], omega_pos[-1]]

power_log = np.log10(dispersion_power)


im = ax.imshow(power_log, extent=extent, origin='lower', aspect='auto',interpolation= "bilinear",cmap='inferno', vmin=-8, vmax=0)

# Scatter dominant points
ax.scatter(dominant_k, dominant_omega, color='cyan', s=50, label='Dominant (k,ω)')

#ax.set_xlim(0.3)
#ax.set_ylim(0,3)
ax.set_xlabel(r"$k$")
ax.set_ylabel(r"$\omega$")
ax.set_title("E(w,k) PSD")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$|E(k,\omega)|^2)$")

plt.tight_layout()
plt.savefig(pjoin(path, "dispersion_psd_positive.pdf"), dpi=300)
plt.show()

