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
DT = metadata_group.attrs['DT']
nE = metadata_group.attrs['nE']
nI = metadata_group.attrs['nI']
nN = metadata_group.attrs['nN']
nB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Tm = metadata_group.attrs['Tm']
Tb = metadata_group.attrs['Tb']
alpha = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mI = metadata_group.attrs['mI']
mN = metadata_group.attrs['mN']
mB = metadata_group.attrs['mB']
density = metadata_group.attrs['density']

# Print the read attributes
print("NC:", NC)
print("NUM_TS:", NUM_TS)
print("write_int:", write_int)
print("write_int_phase:", write_int_phase)
print("DT:", DT)
print("nE:", nE)
print("nI:", nI)
print("nN:", nN)
print("nB:", nB)
print("Te:", Te)
print("Tm:", Tm)
print("Tb:", Tb)
print("alpha:", alpha)
print("beta:", beta)
print("mI:", mI)
print("mN:", mN)
print("mB:", mB)
print("density:", density)
