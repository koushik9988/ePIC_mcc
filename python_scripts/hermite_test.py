import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
from os.path import join as pjoin
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
import matplotlib as mp
from scipy.constants import value as constants


PSI_WEIGHTED = None
DRIFT_VEL = 10  

def init_worker(psi_weighted):
    global PSI_WEIGHTED
    PSI_WEIGHTED = psi_weighted

def compute_hermite_basis(v, M):
    Nv = len(v)
    psi = np.zeros((M, Nv))
    psi[0] = np.pi**(-0.25) * np.exp(-v**2 / 2)
    if M > 1:
        psi[1] = np.sqrt(2) * v * psi[0]
    for m in range(1, M - 1):
        psi[m + 1] = (np.sqrt(2 / (m + 1)) * v * psi[m]- np.sqrt(m / (m + 1)) * psi[m - 1])
    return psi

def get_velocity_range(t, path, species, write_int, drift):
    """Get velocity range AFTER drift subtraction"""
    with h5py.File(path, 'r') as f:
        j = t * write_int
        v = f[f"particle_{species}/{j}"][:, 1] - drift  # Subtract drift here
    return v.min(), v.max()

def process_timestep(t, path, species, write_int, Nx, Nv, v_grid, M_max, n_max, DT_coeff, mfactor, drift):
    global PSI_WEIGHTED
    with h5py.File(path, 'r') as f:
        j = t * write_int
        x = f[f"particle_{species}/{j}"][:, 0]
        v = f[f"particle_{species}/{j}"][:, 1] - drift  # Subtract drift

    x_grid = np.linspace(x.min(), x.max(), Nx)
    ix = np.searchsorted(x_grid, x) - 1
    valid = (ix >= 0) & (ix < Nx)
    ix = ix[valid]
    v = v[valid]

    f_xv = np.zeros((Nx, Nv - 1))
    for i in range(Nx):
        vi = v[ix == i]
        if vi.size:
            f_xv[i], _ = np.histogram(vi, bins=v_grid, density=True)

    f_xv /= Nx

    dx = x_grid[1] - x_grid[0]
    dv = v_grid[1] - v_grid[0]

    # Integral of f(x,v) dxdv should equal 1
    integral_f = np.sum(f_xv) * dx * dv
    if not np.isclose(integral_f, 1.0, atol=1e-3):
        print(f"Warning: Normalization failed at t={t}. Integral = {integral_f}")

    # Hermite projection
    f_m = f_xv @ PSI_WEIGHTED.T

    # Fourier transform
    f_nm = np.fft.fft(f_m, axis=0)[:n_max]

    #f_m_avg = np.mean(f_m, axis=0)
    E_m = np.mean(np.abs(f_m)**2, axis=0)



    time_val = j * DT_coeff * mfactor

    #return t, f_nm, f_m_avg, time_val
    return t, f_nm, E_m, time_val


def main():
    if len(sys.argv) != 3:
        print("Usage: python hermite_fourier_full.py <data_path> <species>")
        sys.exit(1)

    data_path = sys.argv[1]
    species = sys.argv[2]
    path = pjoin(data_path, "result.h5")

    with h5py.File(path, 'r') as f:
        meta = f['/metadata']
        NC = meta.attrs['NC']
        write_int = meta.attrs['write_int_phase']
        NUM_TS = meta.attrs['NUM_TS']
        DT_coeff = meta.attrs['DT_coeff']
        write_interval_phase = meta.attrs['write_int_phase']
        #DT_coeff = meta.attrs['DT_coeff']
        Te = meta.attrs['Te']
        Ti = meta.attrs['Ti']
        alp = meta.attrs['alpha']
        beta = meta.attrs['beta']
        mi = meta.attrs['mI']
        n0 = meta.attrs['density']
        normscheme = meta.attrs.get('norm_scheme', 0)
        EV_TO_K = 11604.52


        eps0 = constants('electric constant')
        kb = constants('Boltzmann constant')
        me = constants('electron mass')
        e = constants('elementary charge')

        # ------------------- Plasma Parameters ---------------------------
        ni0 = n0
        ne0 = n0 / (1 + alp + beta)
        LDi = np.sqrt(eps0 * kb * Ti * EV_TO_K / (ni0 * e**2))
        wpe = np.sqrt(ne0 * e**2 / (eps0 * me))
        wpi = np.sqrt(ni0 * e**2 / (eps0 * mi))

        mfactor = wpi / wpe


    DATA_TS = int(NUM_TS / write_interval_phase) + 1
    ncores = min(cpu_count(), 8)

    # Velocity range scan WITH drift subtraction
    print(f"Scanning velocity range with drift subtraction ({DRIFT_VEL})...")
    with Pool(ncores) as p:
        vr = p.map(partial(get_velocity_range, path=path, species=species, 
                          write_int=write_int, drift=DRIFT_VEL), range(DATA_TS))
    
    vmin = min(v[0] for v in vr)
    vmax = max(v[1] for v in vr)
    
    # Add buffer to ensure all velocities are captured
    v_range = vmax - vmin
    vmin -= 0.1 * v_range
    vmax += 0.1 * v_range
    
    print(f"Velocity range after drift subtraction: [{vmin:.2f}, {vmax:.2f}]")

    Nx = NC
    Nv = 1024
    v_grid = np.linspace(vmin, vmax, Nv)
    v_centers = 0.5 * (v_grid[:-1] + v_grid[1:])
    dv = v_centers[1] - v_centers[0]

    M_max = 250  # hermite modes
    n_max = 250  # k fourier modes

    psi = compute_hermite_basis(v_centers, M_max)
    psi_weighted = psi * dv

    f_nm_time = np.zeros((DATA_TS, n_max, M_max), dtype=complex)
    #f_m_time = np.zeros((DATA_TS, M_max))
    E_m_time = np.zeros((DATA_TS, M_max))


    time_array = np.zeros(DATA_TS)

    worker = partial(process_timestep, path=path, species=species, write_int=write_interval_phase,
                    Nx=Nx, Nv=Nv, v_grid=v_grid, M_max=M_max, n_max=n_max, 
                    DT_coeff=DT_coeff, mfactor=mfactor, drift=DRIFT_VEL)

    print("Processing timesteps...")
    with Pool(ncores, initializer=init_worker, initargs=(psi_weighted,)) as p:
        for t, f_nm, E_m, tt in tqdm(p.imap_unordered(worker, range(DATA_TS)), total=DATA_TS):
            f_nm_time[t] = f_nm
            E_m_time[t] = E_m
            time_array[t] = tt #* write_interval_phase * DT_coeff * mfactor


    #hermite_spec = np.abs(f_m_time)**2
    #Em_x = np.sum(hermite_spec, axis=0)

    hermite_spec = E_m_time

    fourier_spec = np.sum(np.abs(f_nm_time)**2, axis=2)

    outdir = pjoin(data_path, "paper_plots")
    os.makedirs(outdir, exist_ok=True)




    figsize = np.array([150,150/1.618])#Figure size in mm (FOR SINGLE FIGURE)
    dpi = 300                        #Print resolution
    ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

    mp.rc('text', usetex=False)
    mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
    mp.rc('axes', titlesize=10)
    mp.rc('axes', labelsize=10)
    mp.rc('xtick', labelsize=10)
    mp.rc('ytick', labelsize=10)
    mp.rc('legend', fontsize=10)



    # Plot 1: Hermite spectrum over time
    plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)
    plt.contourf(time_array, np.arange(M_max),np.log10(hermite_spec.T + 1e-20), 40, cmap="RdYlBu_r")
    plt.xlabel(r"$\omega_{pi} t$")
    plt.ylabel(r"$m$")
    plt.colorbar(label=r"$\log|f_m|^2$")
    #plt.tight_layout()
    plt.savefig(pjoin(outdir, "fig1_hermite_time.png"), dpi=dpi)
    plt.close()



    plt.figure(figsize=(7, 5))

    #for tidx in [0, DATA_TS//4, DATA_TS//2, DATA_TS-1]:
    for tidx in [0,DATA_TS//6,DATA_TS//4,DATA_TS//3,DATA_TS//2,2*DATA_TS//3,3*DATA_TS//4,DATA_TS-1]:
        plt.plot(np.arange(M_max),hermite_spec[tidx],label=rf"$\omega_{{pi}} t={time_array[tidx]:.1f}$")

    # ----- Reference power-law scaling -----
    tref = DATA_TS //4
    tref1 = DATA_TS//4
    m_range = np.arange(41, M_max) #np.arange(35, M_max) 
    m_range1 = np.arange(24, 43)

    # Normalize by data value at m_range[0]
    f0 = hermite_spec[tref, m_range[0]]
    f1 = hermite_spec[tref1, m_range1[0]]

    #plt.plot(m_range1, f1 * m_range1**(-3/2) / m_range[0]**(-3/2), "r:", label=r"$m^{-3/2}$")
    #plt.plot(m_range1, f1 * m_range1**(-5/2) / m_range1[0]**(-5/2), "k:", label=r"$m^{-5/2}$")
    plt.plot(m_range, f0 * m_range**(-5/2) / m_range[0]**(-5/2), "k:",linewidth=2.5,label=r"$m^{-5/2}$")
    #plt.plot(m_range, f0 * m_range**(-1/2) / m_range[0]**(-1/2), "k--", label=r"$m^{-1/2}$")

    plt.xscale("log")
    plt.yscale("log")
    #remove#plt.ylim(3.5e-11, 2.3e-8)
    plt.xlabel(r"$m$")
    plt.ylabel(r"$\langle |f_m|^2 \rangle_x$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pjoin(outdir, "hermite_power_scaling.pdf"), dpi=dpi)
    plt.show()


    # Plot 3: Fourier-Hermite snapshots
    #snaps = [0, DATA_TS//3, 2*DATA_TS//3, DATA_TS-1]
    snaps = [0,DATA_TS//6,DATA_TS//4,DATA_TS//3,DATA_TS//2,2*DATA_TS//3,3*DATA_TS//4,DATA_TS-1]
    #snaps = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]
    fig, axes = plt.subplots(4, 2, figsize=(10, 8), constrained_layout=True)
    for ax, t in zip(axes.flat, snaps):
        power = np.abs(f_nm_time[t])**2
        #power = np.abs(np.fft.fftshift(f_nm_time[t], axes=0))**2
        im = ax.imshow(np.log10(power.T),origin="lower", aspect="auto", cmap="RdYlBu_r",vmin=-10, vmax=-5)
        #im = ax.imshow(np.log10(np.abs((f_nm_time[t])**2).T),origin="lower", aspect="auto", cmap="RdYlBu_r")
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$m$")
        ax.set_title(rf"$\omega_{{pi}} t = {time_array[t]:.1f}$")
        fig.colorbar(im, ax=ax)
    plt.savefig(pjoin(outdir, "fig5_fourier_hermite.png"), dpi=dpi)
    plt.close()

    
    print(f"Plots saved to {outdir}")

if __name__ == "__main__":
    main()