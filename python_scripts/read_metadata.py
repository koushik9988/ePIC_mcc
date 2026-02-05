import h5py
import os.path
from os.path import join as pjoin
import sys


path = sys.argv[1]
file_name = "result.h5"

file = h5py.File(pjoin(path, file_name), 'r')

metadata_group = file['/metadata']
    
# Read individual attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_int = metadata_group.attrs['write_int']
write_int_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nE = metadata_group.attrs['nE']
nI = metadata_group.attrs['nI']
nN = metadata_group.attrs['nN']
nB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Ti = metadata_group.attrs['Ti']
Tb = metadata_group.attrs['Tb']
alpha = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mI = metadata_group.attrs['mI']
mN = metadata_group.attrs['mN']
mB = metadata_group.attrs['mB']
density = metadata_group.attrs['density']
ve = metadata_group.attrs['ve']
vi = metadata_group.attrs['vi']
vb= metadata_group.attrs['vb']
vn = metadata_group.attrs['vn']
norm_scheme = metadata_group.attrs['norm_scheme']
sub_cycle_interval = metadata_group.attrs['sub_cycle_interval']

# Print the read attributes
print("NC:", NC)
print("NUM_TS:", NUM_TS)
print("write_int:", write_int)
print("write_int_phase:", write_int_phase)
print("DT:", DT_coeff)
print("nE:", nE)
print("nI:", nI)
print("nN:", nN)
print("nB:", nB)
print("Te:", Te)
print("Ti:", Ti)
print("Tb:", Tb)
print("alpha:", alpha)
print("beta:", beta)
print("mI:", mI)
print("mN:", mN)
print("mB:", mB)
print("density:", density)
print("ve:", ve)
print("vi:", vi)
print("vb:", vb)
print("vn:", vn)
print("norm_scheme:",norm_scheme)
print("norm_scheme = 1 : electron_scale \nnorm_scheme = 2 : ion_scale \nnorm_scheme = 3 : subcycling \n")
print("sub_cycle_interval:",sub_cycle_interval)

