
import numpy as np 
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks


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
LDe = metadata.attrs['LDe']
LDi = metadata.attrs['LDi']
density = metadata.attrs['density']
normscheme = metadata.attrs['norm_scheme']

DT = DT_coeff * (1.0 / wpe)

mFactor = wpi/wpe 
if normscheme in [1, 2]:
    mFactor = 1


L = LDe
W = wpe

if normscheme == 2 or normscheme == 4:
    L = LDi
    W = wpi


efield_norm = (density * 1.602176565e-19 * L) / 8.85418782E-12


#====================================================
#  Load Electric Field Data
#====================================================
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for t in time_steps:
    electric_field_data.append(f[f'fielddata/efield/{t}'][:])

EF = np.vstack(electric_field_data)        # (Nt , Nx)
Nt, Nx = EF.shape

x  = np.linspace(0, NC, Nx, endpoint=False)
dx = x[1] - x[0]
time_full = np.arange(Nt) * DT_coeff * write_interval * mFactor


#====================================================
#  Spatial FFT → E(k,t)
#====================================================
Fk = np.fft.fft(EF, axis=1, norm='ortho')
kfft = 2*np.pi*np.fft.fftfreq(NC, dx)

# Only positive k modes
k_modes = kfft[:NC//2]
Fk_pos = Fk[:, :NC//2] * efield_norm * (2/np.sqrt(NC))
energy_spectrum = np.abs(Fk_pos)**2


avg_power = np.mean(energy_spectrum, axis=0)

peaks, _ = find_peaks(avg_power, height=np.max(avg_power)*0.001, distance=2)

print("Detected modes at k:", k_modes[peaks])

dt_snapshot = DT_coeff * write_interval * mFactor
phase_labels = []

for idx in peaks:
    k_dim = k_modes[idx]              
   
    Ek_t = Fk_pos[:, idx]                

    wspec = np.fft.fft(Ek_t)
    freqs = 2*np.pi*np.fft.fftfreq(Nt, dt_snapshot)

    # keep positive ω only
    pos = slice(0, Nt//2)
    peak_ω = freqs[pos][ np.argmax( np.abs(wspec[pos]) ) ]

    vph = peak_ω/k_dim                  

    #phase_labels.append(f"v_ph ≈ {vph:.3e}")
    phase_labels.append(f"$v_\phi$ = {vph:.3f}")



plt.figure()

for j, idx in enumerate(peaks):
    Ek = energy_spectrum[:, idx]
    plt.plot(time_full, np.log10(Ek), label=phase_labels[j])

plt.xlabel(r"$\omega_{pi}t$")
#plt.ylim([np.max(np.log10(energy_spectrum))*0.6, np.max(np.log10(energy_spectrum))*1.1])
plt.ylabel(r"$|E_k|^2$")
#plt.title("Evolution of Spectral Modes vs Time\n(labelled by Phase Velocity)")
plt.grid(alpha=0.4)
plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig(pjoin(path,"Ek_vs_time_phase_velocity.pdf"),dpi=300)
plt.show()
