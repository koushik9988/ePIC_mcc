import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py
import matplotlib.animation as animation


#hdf5 file name and path 
file_name = 'result.h5'
path = sys.argv[1]

plot_path = './plots'

path_fig = pjoin(path,plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
#------------------------------------------

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = e*metadata_group.attrs['Te']
Ti = e*metadata_group.attrs['Ti']
Tb = e*metadata_group.attrs['Tb']
Tn = e*metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS /write_interval_phase) + 1

#-----------------------Debye lenght calculation--------
ni0 = n0
ne0 = n0*(1- alp - alp*beta)
ni0 = n0
nn0 = alp*ni0  
nb0 = beta*ni0

LD = np.sqrt(eps0*Te/(ne0*e**2)) 
#---------------------------------------------------------

fig, (ax_pot, ax_phase1, ax_phase2, ax_phase3, ax_phase4) = plt.subplots(5, 1, figsize=(10, 8))

def animate(i):
    j = i*write_interval_phase

    data_phase_e = f["particle_electron/%d"%j]
    
    dataex = data_phase_e[:,0]
    dataevx = data_phase_e[:,1]


    data_phase_i = f["particle_ion/%d"%j]
    
    dataix = data_phase_i[:,0]
    dataivx = data_phase_i[:,1]

    data_phase_n = f["particle_negion/%d"%j]
    
    datanx = data_phase_n[:,0]
    datanvx = data_phase_n[:,1]

    data_phase_b = f["particle_beam/%d"%j]
    
    databx = data_phase_b[:,0]
    databvx = data_phase_b[:,1]


    ax_phase1.clear()
    ax_phase1.scatter(dataex,dataevx,marker='.',color='b',alpha=1.0,s=10, label = "electron phase space")
    ax_phase1.set_xlabel('$x$')
    ax_phase1.set_ylabel('$v$')
    ax_phase1.legend(loc='upper right',framealpha=0.5)

    ax_phase2.cla()
    ax_phase2.scatter(dataix,dataivx,marker='.',color='r',alpha=1.0,s=10, label = "ion phase space")
    ax_phase2.set_xlabel('$x$')
    ax_phase2.set_ylabel('$v$')
    ax_phase2.legend(loc='upper right',framealpha=0.5)

    ax_phase3.cla()
    ax_phase3.scatter(datanx,datanvx,marker='.',color='g',alpha=1.0,s=10, label = "negative-ion phase space")
    ax_phase3.set_xlabel('$x$')
    ax_phase3.set_ylabel('$v$')
    ax_phase3.legend(loc='upper right',framealpha=0.5)

    ax_phase4.cla()
    ax_phase4.scatter(databx,databvx,marker='.',color='y',alpha=1.0,s=10, label = "beam phase space")
    ax_phase4.set_xlabel('$x$')
    ax_phase4.set_ylabel('$v$')
    ax_phase4.legend(loc='upper right',framealpha=0.5)

    ax_pot.clear()
    pot = f["fielddata/pot/%d"%j]
    #print(pot[:])
    x = np.linspace(0,NC, len(pot))
    ax_pot.plot(x, pot[:], color='black', label="Potential")
    ax_pot.set_ylabel('$\phi$')
    ax_pot.legend(loc='upper right',framealpha=0.5)

    if(save_fig == 1):
        if(i%1000 == 0):
            plt.savefig(pjoin(path_fig, 'phase-pot_%d.png' % (j)), dpi=1200)

    return ax_phase1, ax_phase2 , ax_phase3, ax_phase4 ,ax_pot


ani = animation.FuncAnimation(fig, animate, frames = DATA_TS_PHASE, blit= False, interval= 100, repeat=False)

plt.show()
