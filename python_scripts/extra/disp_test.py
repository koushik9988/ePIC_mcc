import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

# Input file paths
file_name = 'result.h5'

path1 = sys.argv[1]
path2 = sys.argv[2]


f = h5py.File(pjoin(path1, file_name), 'r')
    
    #----------------Read metadata ------------
metadata_group = f['/metadata']
metadata_electron = f['/metadata_species/electron']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']
NUM_TS = metadata_group.attrs['NUM_TS']
DT_coeff = metadata_group.attrs['DT_coeff']
write_interval  = metadata_group.attrs['write_int']
DT = DT_coeff * (1.0 / we)
NC = metadata_group.attrs['NC']
dx = 1

wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)

NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*LD) # Normalized 

def load_data(path, file_name):
    f = h5py.File(pjoin(path, file_name), 'r')
    
    #----------------Read metadata ------------
    metadata_group = f['/metadata']
    metadata_electron = f['/metadata_species/electron']

    # Read necessary constants and metadata
    eps0 = constants('electric constant')
    kb = constants('Boltzmann constant')
    me = constants('electron mass')
    AMU = constants('atomic mass constant')
    e = constants('elementary charge')

    NC = metadata_group.attrs['NC']
    NUM_TS = metadata_group.attrs['NUM_TS']
    write_interval  = metadata_group.attrs['write_int']
    DT_coeff = metadata_group.attrs['DT_coeff']
    LD = metadata_group.attrs['LDe']
    LDi = metadata_group.attrs['LDi']
    we = metadata_group.attrs['wpe']
    wp = metadata_group.attrs['wpi']
    n0 = metadata_group.attrs['density']
    Te = metadata_electron.attrs['temperature']
    mi = metadata_electron.attrs['mass']

    EV_TO_K = 11604.52 
    DT = DT_coeff * (1.0 / we)

    # Read electric field data
    electric_field_data = []
    time_steps = sorted(map(int, f['fielddata/efield'].keys()))
    for time_step in time_steps:
        EF_data = f['fielddata/efield/' + str(time_step)]
        electric_field_data.append(EF_data[:])  # Append electric field data

    # Combine electric field data into a 2D array
    EF = np.vstack(electric_field_data)
    x = np.linspace(0, NC, EF.shape[1])
    dx = x[1] - x[0]

    # Reduce EF data based on time steps
    wpet_1 = 0
    wpet_2 = NUM_TS * DT_coeff
    y1 = wpet_1 / (DT_coeff * write_interval)
    y2 = wpet_2 / (DT_coeff * write_interval)
    E = EF[int(y1):int(y2), :]

    # Compute FFT
    F = np.fft.fftn(E, norm='ortho')

    NUM_TS1 = wpet_1/DT_coeff
    NUM_TS2 = wpet_2/DT_coeff
    actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT)
    omega = 2 * np.pi * np.arange(NUM_TS) / actual_sim_time  # Unnormalized Omega
    k = 2 * np.pi * np.arange(NC + 1) / (NC * dx * LD)  # Normalized 

    # Meshgrid for Omega and K
    Omega, K = np.meshgrid(omega, k, indexing='ij')

    # Halve data (since it's symmetric in FFT)
    halflen = np.array(F.shape, dtype=int) // 2
    Omega = Omega[:halflen[0], :halflen[1]]
    K = K[:halflen[0], :halflen[1]]
    F = F[:halflen[0], :halflen[1]]

    # Normalize Omega and K
    Omega /= we  # Normalize Omega
    K = K * LD  # Normalize K

    Z = np.log(np.abs(F))
    return K, Omega, Z

# Load data for both paths (datasets)
K1, Omega1, Z1 = load_data(path1, file_name)
K2, Omega2, Z2 = load_data(path2, file_name)


Z_diff = np.abs(Z1 - Z2)



EV_TO_K = 11604.52
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

Ti = 0.026
Te = 1
mi = AMU

alpha1 = 1/(1+1+0.4)
raw_analytic_ILW1 = True
if raw_analytic_ILW1:    
    ilw1 = np.sqrt(((Ti*EV_TO_K*kb/mi)*(1 + alpha1 + 3*(k**2*LDi**2) + 3*(1-alpha1)*(Ti/Te)))/((k*LDi)**2 + (Ti/Te)*(1-alpha1)))*k

alpha2 = 20/(1+20+0.4)
raw_analytic_ILW2 = False
if raw_analytic_ILW2:    
    ilw2 = np.sqrt(((Ti*EV_TO_K*kb/mi)*(1 + alpha2 + 3*(k**2*LDi**2) + 3*(1-alpha2)*(Ti/Te)))/((k*LDi)**2 + (Ti/Te)*(1-alpha2)))*k

# Plotting parameters
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))

# Plot dispersion from both datasets on the same axes (first subplot)
c1 = ax1.pcolor(K1, Omega1, Z1, cmap='rainbow', shading='auto', vmin= -10,vmax=10)
#c1 = ax1.contourf(K1, Omega1, Z1, cmap='rainbow',shading='auto', vmin=np.min(Z1),vmax=np.max(Z1))
c2 = ax1.pcolor(K2, Omega2, Z2, cmap='rainbow', shading='auto', vmin=-10,vmax=10)
#c2 = ax1.contourf(K2, Omega2, Z2, cmap='rainbow',shading='auto',vmin=np.min(Z1),vmax=np.max(Z1))

if raw_analytic_ILW1:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    ax1.plot(k*LD, ilw1/we, color='k', linestyle='-.', lw = 1.0, label='$IAW1$')
if raw_analytic_ILW2:
    #plt.plot(k*LD, ep, color='k', linestyle='--', lw = 1.0, label='$EPW$')
    ax1.plot(k*LD, ilw2/we, color='k', linestyle='-.', lw = 1.0, label='$IAW2$')

# Colorbar for the first plot
cbar1 = plt.colorbar(c2, ax=ax1)
cbar1.set_label('$\zeta$')

# Labels and limits for the first plot
ax1.set_xlabel('$k \lambda_{D}$')
ax1.set_ylabel('$\omega/\omega_{pe}$')
ax1.set_xlim([0, 3.5])
ax1.set_ylim([0, 1])
ax1.set_title('Dispersion Comparison')

# Plot the absolute difference in a separate plot (second subplot)
c3 = ax2.pcolor(K1, Omega1, Z_diff, cmap='viridis', shading='auto')

# Colorbar for the difference plot
cbar2 = plt.colorbar(c3, ax=ax2)
cbar2.set_label('Absolute Difference |$\zeta_1 - \zeta_2$|')

# Labels and limits for the second plot
ax2.set_xlabel('$k \lambda_{D}$')
ax2.set_ylabel('$\omega/\omega_{pe}$')
ax2.set_xlim([0, 3.5])
ax2.set_ylim([0, 2])
ax2.set_title('Absolute Difference in Dispersion')

# Save and show the combined plot
plt.tight_layout()
plt.savefig(pjoin(path1, 'dispersion_comparison_with_difference.png'), dpi=1200)
plt.show()
