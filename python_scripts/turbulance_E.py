import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq
from scipy.signal import correlate
import matplotlib as mp

if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_data>")
    sys.exit(1)

path = sys.argv[1]

output_dir = './turbulance_plots'
path_fig = pjoin(path, output_dir)
os.makedirs(path_fig, exist_ok=True)

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

mfactor = 1#wpi/wpe

# Load Electric Field Data 
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
for time_step in time_steps:
    EF_data = f[f'fielddata/efield/{time_step}']
    electric_field_data.append(EF_data[:])

EF = np.vstack(electric_field_data)  # Shape: (Nt, Nx)
x = np.linspace(0, NC, EF.shape[1])
dx = x[1] - x[0]

print(EF.shape[0])

time_full = np.arange(EF.shape[0]) * DT_coeff * write_interval * mfactor #normalized
#(np.arrange(n)=> 0,1,2,3....n-1)

#1:Compute fluctuations/
mean_E_over_x = np.mean(EF, axis=1)[:, np.newaxis]  # Mean over space for each time: shape (Nt, 1)
fluct_E = EF #- mean_E_over_x 

#mean_E_over_t = np.mean(EF, axis=0)[np.newaxis, :]  # Mean over time for each position: (1, Nx)
#fluct_E = EF - mean_E_over_t

figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)


plt.figure(figsize = figsize/10.4,constrained_layout=True,dpi=ppi)
plt.imshow(EF, aspect='auto', cmap='viridis', origin='lower', extent=[0, NC, 0, time_full[-1]])
plt.xlabel('$x$')
plt.ylabel('$\omega_{pi}t$')
plt.title('$E(x,t)$')
plt.colorbar(label='$E(x,t)$')
plt.savefig(pjoin(path_fig, 'efield_contour.png'))
plt.close()

#Spatial autocorrelation at a fixed time step
t_step = 4000 #max is 3001 for our dataset #fix time later
#correlate(in1, in2, mode='full', method='auto')
#in1 =>array_like First input.
#in2 =>array_like Second input. Should have the same number of dimensions as in1.

#modestr {‘full’, ‘valid’, ‘same’}, optional
#A string indicating the size of the output:

#full => The output is the full discrete linear cross-correlation of the inputs. (Default)

#valid=>The output consists only of those elements that do not rely on the zero-padding. In ‘valid’ mode, either in1 or in2 must be at least as large as the other in every dimension.

#same=>The output is the same size as in1, centered with respect to the ‘full’ output.


auto_corr = correlate(fluct_E[t_step, :], fluct_E[t_step, :], mode='full')
auto_corr = auto_corr[auto_corr.size//2:] / auto_corr[auto_corr.size//2]  # Normalize and take positive lags (lags= amount of shift applied to the signal. here it of 1 debye lenght(normalized so 1))
lags = np.arange(len(auto_corr))*dx

plt.figure(figsize = figsize/10.4,constrained_layout=True,dpi=ppi)
plt.plot(lags, auto_corr)
plt.xlabel('$(Δx)$')
plt.ylabel(r'$\sum E(x)E(x+\Delta x)$')
plt.title(f'Spatial Autocorrelation at time  ={t_step*DT_coeff*write_interval*mfactor:.2f}')
plt.grid(True)
plt.savefig(pjoin(path_fig, 'spatial_autocorrelation.png'))
plt.close()

# Temporal autocorrelation at fixed spatial loacation
x_pos = NC-1
auto_corr_temp = correlate(fluct_E[:, x_pos], fluct_E[:, x_pos], mode='full')
auto_corr_temp = auto_corr_temp[auto_corr_temp.size//2:] / auto_corr_temp[auto_corr_temp.size//2]
lags_temp = np.arange(len(auto_corr_temp)) * write_interval * DT_coeff *  mfactor 

plt.plot(lags_temp, auto_corr_temp)
plt.xlabel('$Δt$')
plt.ylabel(r'$\sum E(t)E(t + \Delta t )$')
plt.title(f'Temporal Autocorrelation at x={x_pos * dx:.2f}')
plt.savefig(pjoin(path_fig, 'temporal_autocorrelation.png'))

#print("len(auto_corr_temp))",len(auto_corr_temp))
# Spatial power spectrum at a fixed time step
#(time averaging also works)
Nx = EF.shape[1]
fft_E = fft(fluct_E[t_step, :])
k_freq = fftfreq(Nx, d=dx)[:Nx//2]  #Positive wave numbers
power_spectrum = np.abs(fft_E[:Nx//2])**2

#Fit power law (kolmogorov like)
def power_law(k, alpha, A):
    return A * k**alpha

#popt, pcov = curve_fit(power_law, k_freq, power_spectrum)
#popt → an array of the best-fit parameters
#pcov → the covariance matrix of the parameters

fit_start = 20
popt, _ = curve_fit(power_law, k_freq[fit_start:], power_spectrum[fit_start:])
print(popt) 
#curve fit gives two sets of value we only need fit parameters
print(f'Fitted power law slope (alpha): {popt[0]:.2f}')

plt.figure(figsize = figsize/10.4,constrained_layout=True,dpi=ppi)
plt.loglog(k_freq[1:], power_spectrum[1:])
plt.loglog(k_freq[fit_start:], power_law(k_freq[fit_start:], *popt), 'r--', label=f'Fit: k^{popt[0]:.2f}')
plt.xlabel('$k$')
plt.ylabel('$|E(k)|^2$')
plt.legend()
plt.grid(True)
plt.savefig(pjoin(path_fig, 'power_spectrum_fit.png'))
plt.close()

print("end...")