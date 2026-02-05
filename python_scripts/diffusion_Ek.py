"""
Simple python script to analyze Electric field to study turbulence.
@kaushik kalita
@Date : 16/10/2025
Works with datasets of 1st paper/work.
"""
import numpy as np
import h5py
import sys
from os.path import join as pjoin
from scipy.constants import value as constants 
from scipy.stats import linregress  
import os
from os.path import join as pjoin
import matplotlib.pyplot as plt


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

#----------------------------------------------

#-------------------------norm------------------
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

mfactor = wpi/wpe
if normscheme == 2 or normscheme == 1:
    mfactor = 1
print("mfactor:",mfactor)

DT = DT_coeff * (1.0 / wpi)

efield_norm = (density * 1.602176565e-19 * LDi) / 8.85418782E-12
#----------------------------------------------

#----------------Load Electric Field Data---------------------------------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)

###

print(f"Sliced Shape (Turbulent Region): {EF.shape}")


Nt, Nx = EF.shape
x = np.linspace(0, NC, Nx)
#normalized unit
dx = x[1] - x[0]
#----------------------------------------------


#normalized k
kfft = 2 * np.pi * np.fft.fftfreq(Nx, dx)
k_modes = kfft[:Nx // 2]
F_k = np.fft.fft(EF, axis=1, norm='ortho')[:, :Nx // 2]


power_spectrum = np.mean(np.abs(F_k)**2, axis=0)
max_idx = np.argmax(power_spectrum)
k_val = abs(k_modes[max_idx])

#un-normalize
k_val  = k_val/LDi

print(f"Detected most unstable mode index: {max_idx}")
print("k_max : ",k_val*LDi)

#un-normalize the field
Ek = np.abs(F_k[:, max_idx])*efield_norm*(2/np.sqrt(NC)) #scaled to match physical value( extra factor we take -k modes ??)

"""
fft result are proportional to have to scale by 1/N or in case of ortho(where persaval idenity is preserved) 1/sqrt(N) 
and a multiplication factor of 2 comes because energy is divided into both positive and negative spectra , so positive spectra carry half and negative half.

"""

Ek_mean = np.mean(Ek)
Ek_max = np.max(Ek)
Ek2_max = np.max(Ek**2)
Ek2_mean = np.mean(Ek**2)

print("Max field amplitude in volts",Ek_max)

#Resonace velocity widthmmm
delta_vk_half = 2.0 * np.sqrt((e * Ek2_max) / (species_mass * k_val))
delta_vk = 2*delta_vk_half

print("resonance velocity width",delta_vk/(LDi*wpi))


psdk = np.mean(np.abs(Ek)**2, axis=0) #time averaged

# (pi*e^2/2*m*2)*sum(E_k^2 delta(omega - k v))
Dth = (np.pi * e**2 * psdk) / (2 * species_mass **2 * k_val * delta_vk)
#Dth =  (e**2 * psdk) / (mb**2 * k_val * delta_vk)
#**********************************************************************************************************


#*************************velocity variance***************************************************************
data = f["time_var/kinetic_energy"]
ts = data[:, 0] * mfactor  # ω_pi * t

# Reference velocity: initial mean (t=0)
data_phase0 = f[f"particle_{particle_type}/0"]
v0 = data_phase0[:, 1]*(LDi*wpi)
v0_mean = np.mean(v0)
print("v_mean =", v0_mean)

#empty lists to store results
sigma_v2 = []
t_vals = []

for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    data_phase_t = f[f"particle_{particle_type}/{j}"]
    v_t = data_phase_t[:, 1]*(LDi*wpi)
    sigma2 = np.mean((v_t - v0_mean) ** 2)
    sigma_v2.append(sigma2)
    t_vals.append(j*DT*mfactor)

sigma_v2 = np.array(sigma_v2)
t_vals = np.array(t_vals)

#*********************************************************************************************************

#staring point for fit
t_start = 3.7e-5

#Filter out all points before t_start
mask = t_vals >= t_start
t_fit = t_vals[mask]
sigma_fit = sigma_v2[mask]


coeff = np.polyfit(t_fit, sigma_fit, 1)
fit_line = np.polyval(coeff, t_fit)

D_fit = coeff/2


#*****************
print("difffusion coeffiecnet",D_fit[0])
print("diffusion theory E_k",Dth)
#****************

# Plot
plt.plot(t_vals, sigma_v2, label=r'$\sigma_v^2(t)$')
plt.plot(t_fit, fit_line, 'r--', label=f'slope={coeff[0]:.3e}')

plt.legend()
#plt.xlabel(r'$\omega_{pi} t$')
plt.xlabel(r'$t$')
plt.ylabel(r'$\sigma_{v}^2$')

textstr = '\n'.join((r'$D_{\mathrm{theory}} = %.3e$' % Dth, r'$D_{\mathrm{fit}} = %.3e$' % D_fit[0]))
plt.text(0.65, 0.75, textstr, transform=plt.gca().transAxes,fontsize=11, verticalalignment='top',bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.5))

plt.savefig(pjoin(path_fig, f"velocity_diffusion_{particle_type}.png"))

plt.show()

