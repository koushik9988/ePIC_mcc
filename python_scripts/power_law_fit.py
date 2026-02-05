import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.optimize import curve_fit

if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]

file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

# Read Metadata 
metadata = f['/metadata']
NC = metadata.attrs['NC']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
wpi = metadata.attrs['wpi']


DT = DT_coeff * (1.0 / wpe)

# Load Electric Field Data 
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    #EF_data = f[f'fielddata/efield/{time_step}']
    EF_data = f[f'fielddata/den_negion/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
x = np.linspace(0, NC, EF.shape[1])
dx = x[1] - x[0]

print("Shape of EF:", EF.shape)
print("(x,t) = ",EF.shape[1],EF.shape[0])

time_full = np.arange(EF.shape[0]) * DT_coeff * write_interval

# Spatial FFT
F_k = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)
k_modes = kfft[:NC // 2]    #Positive k modes
F_k_pos = F_k[:, :NC // 2]  # Positive k modes only

# Power Spectrum mean for all time steps
power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)

# curve fit range
k_min = float(input("Enter minimum k for fitting: "))
k_max = float(input("Enter maximum k for fitting: "))
#check
fit_range = (k_modes >= k_min) & (k_modes <= k_max)
k_fit = k_modes[fit_range]
P_fit = power_spectrum[fit_range]

#power law turbulance
def power_law(k, A, alpha):
    return A * k**(-alpha)

#logx and logy fit
log_k = np.log(k_fit)
log_P = np.log(P_fit)
popt, pcov = curve_fit(lambda k, alpha, logA: alpha * k + logA, log_k, np.log(P_fit))
alpha_fit = -popt[0]  # slope = -alpha
A_fit = np.exp(popt[1])

print(f"Fitted power law: |E(k)|^2 ~ k^-{alpha_fit:.3f}")

fig, ax = plt.subplots(figsize=(8, 5))
ax.loglog(k_modes, power_spectrum, label=r'$|E(k)|^2$')
ax.loglog(k_fit, power_law(k_fit, A_fit, alpha_fit), 'r--',label=rf'Fit: $k^{{-{alpha_fit:.2f}}}$')
ax.set_xlabel('$k$')
ax.set_ylabel('$|E(k)|^2$')
ax.grid(True, which='both', ls='--', lw=0.5)
ax.legend()
plt.tight_layout()
plt.savefig(pjoin(path, 'turb.png'), dpi=300)
plt.show()
