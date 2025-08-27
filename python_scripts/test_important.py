import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin

# --------------------- Argument Handling ---------------------
if len(sys.argv) != 3:
    print("Usage: python plot_power_spectrum_at_time.py <path_to_data> <time_step_index>")
    sys.exit(1)

path = sys.argv[1]
time_index = int(sys.argv[2])

# --------------------- Load File and Metadata ---------------------
file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

metadata = f['/metadata']
NC = metadata.attrs['NC']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']

DT = DT_coeff * (1.0 / wpe)

# --------------------- Get List of Time Steps ---------------------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
if time_index < 0 or time_index >= len(time_steps):
    print(f"Error: time_index {time_index} out of bounds (0 to {len(time_steps)-1})")
    sys.exit(1)

ts = time_index*write_interval
EF = f[f'fielddata/efield/{ts}'][:]  # Shape: (NC,)

x = np.linspace(0, NC, EF.shape[0], endpoint=False)
dx = x[1] - x[0]

# --------------------- FFT at Selected Time ---------------------
F_k = np.fft.fft(EF, norm='ortho')
kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)
k_modes = kfft[:NC // 2]
power_spectrum = np.abs(F_k[:NC // 2])**2

# --------------------- Plot ---------------------
plt.figure(figsize=(8, 5))
plt.plot(k_modes, power_spectrum, label=f'Time step index = {time_index} (t = {ts:.2f})')
plt.xlabel('$k$')
plt.ylabel('$|E(k)|^2$')
plt.grid(True)
plt.legend()
plt.tight_layout()

output_name = f'power_spectrum_t{time_index}.pdf'
plt.savefig(pjoin(path, output_name), dpi=300)
plt.show()
