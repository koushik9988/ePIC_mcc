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


PSI_WEIGHTED = None

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

def get_velocity_range(t, path, species, write_int):
    with h5py.File(path, 'r') as f:
        j = t * write_int
        v = f[f"particle_{species}/{j}"][:, 1]
    return v.min(), v.max()

def process_timestep(t, path, species, write_int,Nx, Nv, v_grid,M_max, n_max, DT_coeff, mfactor):
    global PSI_WEIGHTED
    with h5py.File(path, 'r') as f:
        j = t * write_int
        x = f[f"particle_{species}/{j}"][:, 0]
        v = f[f"particle_{species}/{j}"][:, 1]

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

    # Hermite projection
    f_m = f_xv @ PSI_WEIGHTED.T

    # Fourier transform
    f_nm = np.fft.fft(f_m, axis=0)[:n_max]

    f_m_avg = np.mean(f_m, axis=0)
    time_val = j * DT_coeff * mfactor

    return t, f_nm, f_m_avg, time_val

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

    DATA_TS = int(NUM_TS / write_int) + 1
    ncores = min(cpu_count(), 16)

    # Velocity range scan
    with Pool(ncores) as p:
        vr = p.map(partial(get_velocity_range,path=path,species=species,write_int=write_int),range(DATA_TS))
    vmin = min(v[0] for v in vr)
    vmax = max(v[1] for v in vr)
    vmin -= 0.5 * abs(vmin)
    vmax += 0.5 * abs(vmax)

    Nx = NC
    Nv = 1024
    v_grid = np.linspace(vmin, vmax, Nv)
    v_centers = 0.5 * (v_grid[:-1] + v_grid[1:])
    dv = v_centers[1] - v_centers[0]

    M_max = 500 # hermite modes , hermite polynomial orders basically
    n_max = 100 # k fourier modes

    psi = compute_hermite_basis(v_centers, M_max)
    psi_weighted = psi * dv

    f_nm_time = np.zeros((DATA_TS, n_max, M_max), dtype=complex)
    f_m_time = np.zeros((DATA_TS, M_max))
    time_array = np.zeros(DATA_TS)

    worker = partial(process_timestep,path=path,species=species,write_int=write_int,Nx=Nx,Nv=Nv,v_grid=v_grid,M_max=M_max,n_max=n_max,DT_coeff=DT_coeff,mfactor=1.0)

    with Pool(ncores,initializer=init_worker,
              initargs=(psi_weighted,)) as p:
        for t, f_nm, f_m, tt in tqdm(p.imap_unordered(worker, range(DATA_TS)),total=DATA_TS):
            f_nm_time[t] = f_nm
            f_m_time[t] = f_m
            time_array[t] = tt

    hermite_spec = np.abs(f_m_time)**2
    fourier_spec = np.sum(np.abs(f_nm_time)**2, axis=2)

    outdir = pjoin(data_path, "paper_plots")
    os.makedirs(outdir, exist_ok=True)

    
    plt.figure(figsize=(7,5))
    plt.contourf(time_array, np.arange(M_max),
                 np.log10(hermite_spec.T + 1e-20), 40,
                 cmap="RdYlBu_r")
    plt.xlabel(r"$\omega t$")
    plt.ylabel("Hermite mode $m$")
    plt.colorbar(label=r"$\log_{10}|f_m|^2$")
    plt.tight_layout()
    plt.savefig(pjoin(outdir, "fig3_hermite_time.png"), dpi=300)
    #plt.close()

    
    plt.figure(figsize=(7,5))
    for tidx in [DATA_TS//4, DATA_TS//2, DATA_TS-1]:
        plt.plot(hermite_spec[tidx], label=rf"$t={time_array[tidx]:.1f}$")
    m = np.arange(1, M_max)
    plt.plot(m, m**(-3/2), "k--", label=r"$m^{-3/2}$")
    plt.plot(m, m**(-5/2), "k:", label=r"$m^{-5/2}$")
    plt.plot(m, m**(-1/2), "k:", label=r"$m^{-1/2}$")
    plt.yscale("log")
    plt.xlabel("Hermite mode $m$")
    plt.ylabel(r"$\langle |f_m|^2 \rangle$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pjoin(outdir, "fig4_hermite_cuts.png"), dpi=300)
    #plt.close()

    
    snaps = [0, DATA_TS//3, 2*DATA_TS//3, DATA_TS-1]
    fig, axes = plt.subplots(2,2,figsize=(10,8), constrained_layout=True)
    for ax, t in zip(axes.flat, snaps):
        im = ax.imshow(np.log10(np.abs(f_nm_time[t])**2 + 1e-20).T,
                       origin="lower", aspect="auto", cmap="RdYlBu_r")
        ax.set_xlabel("Fourier mode $n$")
        ax.set_ylabel("Hermite mode $m$")
        ax.set_title(rf"$t={time_array[t]:.1f}$")
        fig.colorbar(im, ax=ax)
    plt.savefig(pjoin(outdir, "fig5_fourier_hermite.png"), dpi=300)
    #plt.close()

    
    F = np.sum(hermite_spec[:,1:], axis=1)
    plt.figure(figsize=(7,5))
    plt.plot(time_array, F)
    plt.xlabel(r"$\omega t$")
    plt.ylabel("Free energy (Hermite)")
    plt.tight_layout()
    plt.savefig(pjoin(outdir, "fig7_free_energy.png"), dpi=300)
    #plt.close()

    plt.show()

if __name__ == "__main__":
    main()
