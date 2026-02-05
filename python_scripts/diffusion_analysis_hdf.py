import h5py
import numpy as np
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys, os
import matplotlib as mp


if len(sys.argv) != 3:
    print("Usage: python analyze_all_diffusion_plots.py <path_to_simulation>")
    sys.exit(1)


path = sys.argv[1]
ptype = sys.argv[2]

analysis_dir = pjoin(path, "new_analysis")
fname = pjoin(analysis_dir, f"local_velocity_diffusion_{ptype}.h5")


path_fig = pjoin(path, "diffusion_analysis_plots")
os.makedirs(path_fig, exist_ok=True)

with h5py.File(fname, "r") as f:

    axes = f["axes"]
    meta = f["metadata"].attrs

    v_centers = axes["v_centers"][:]
    tau_phys  = axes["tau_phys"][:]
    tau_steps      = axes["tau_steps"][:]  

    all_pdfs = f["pdfs/P_dv_tau_v"][:]           # [vbin, tau, dv]
    sigma    = f["gaussian_fit/sigma"][:]        # [vbin, tau]
    mu       = f["gaussian_fit/mu"][:]

    D_var   = f["local_diffusion/D_from_var"][:]
    D_msd   = f["local_diffusion/D_from_msd"][:]
    D_sigma = f["local_diffusion/D_from_sigma"][:]

    pdf_centers_grp = axes["pdf_centers"]
    pdf_centers = {k: pdf_centers_grp[k][:] for k in pdf_centers_grp.keys()}

    v_phase = meta["v_phase"]
    Nfit    = meta["Nfit"]

    all_var = f["statistics/variance"][:]
    all_sigma_fit = f["gaussian_fit/sigma"][:] # Gaussian σ


nv_bins, n_tau, _ = all_pdfs.shape


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






"""
#data_run_16
vphi = np.array([
 0.109198,
 0.185772,
 0.200713,
 0.203467,
 0.206778,
 0.215871,
 0.221355,
 0.315984,
 0.341824,
 0.367474,
 0.367777,
 0.368579,
 0.460106,
 0.48307,
 0.544051,
 0.546292,
 0.563223,
 0.57875,
 0.608893,
 0.634263,
 0.648226,
 0.655189,
 0.683705,
 0.714579,
 0.761279,
 0.923303,
 0.938633,
 0.957222,
 1.07444,
 1.09291,
 1.15038,
 1.16799,
 1.217,
 1.29384,
 1.33779,
 1.389,
 1.61736,
 1.62989,
 1.65882,
 1.7411,
 1.75316,
 1.7612,
 1.78097,
 1.79177,
 1.82306,
 1.83594,
 1.91249,
 1.91368,
 1.98861,
 2.0452,
 2.16337,
 2.27971,
 2.56425,
 2.71375,
 2.73977,
 2.79856,
 2.89504,
 2.89801,
 2.92788,
 3.03844,
 3.11694,
 3.20696,
 3.20914,
 3.33653,
 3.41282,
 3.53616,
 3.63527,
 3.70814,
 3.79167,
 3.88683,
 3.9676,
 4.13427,
 4.30321,
 4.43364,
 4.64365,
 4.89623,
 5.03694,
 5.48497,
 6.1214,
 6.41709,
 6.74657,
 7.116,
 8.06758,
 8.91498,
 10.0589,
 11.0169,
 15.279,
 18.9541,
 22.5713,
 42.8854
])

Dth = np.array([
 0.000431776,
 0.00025826,
 0.000241214,
 0.000290633,
 0.000233728,
 0.000346016,
 0.000239247,
 0.000269867,
 0.000253761,
 0.000417268,
 0.000364332,
 0.000285624,
 0.00029888,
 0.000263697,
 0.00029791,
 0.000225924,
 0.000318127,
 0.000229526,
 0.000289504,
 0.000240168,
 0.000283167,
 0.000304115,
 0.000397314,
 0.000442979,
 0.000347444,
 0.00024861,
 0.000405453,
 0.000375195,
 0.000409805,
 0.000410569,
 0.000393169,
 0.000301447,
 0.000323618,
 0.000464232,
 0.000215771,
 0.000459083,
 0.000266326,
 0.000509729,
 0.00037482,
 0.000394469,
 0.00045922,
 0.000417466,
 0.000406118,
 0.00044567,
 0.000562002,
 0.000289246,
 0.00051048,
 0.000508298,
 0.000532651,
 0.00049126,
 0.00051508,
 0.00049665,
 0.000595442,
 0.000584763,
 0.000565665,
 0.000544212,
 0.000611426,
 0.00067192,
 0.000572181,
 0.000658104,
 0.000832026,
 0.000772701,
 0.000805294,
 0.000916015,
 0.000988612,
 0.00099738,
 0.00102893,
 0.00138147,
 0.0011757,
 0.0015322,
 0.00154003,
 0.00177463,
 0.00227457,
 0.00162843,
 0.00223317,
 0.0021027,
 0.00215151,
 0.00233751,
 0.00315816,
 0.00286015,
 0.00391872,
 0.0037279,
 0.00370321,
 0.00481947,
 0.0047832,
 0.00428303,
 0.00509821,
 0.00755132,
 0.00863812,
 0.031345
])
"""



""""
#alpha = 0
vphi = np.array([
 0.0208706,
 0.0543248,
 0.14269,
 0.199083,
 0.208796,
 0.21415,
 0.222759,
 0.232089,
 0.343256,
 0.360399,
 0.374085,
 0.381455,
 0.402673,
 0.476657,
 0.498409,
 0.503739,
 0.527727,
 0.581303,
 0.634922,
 0.718295,
 0.732449,
 0.777703,
 0.808082,
 0.820909,
 0.838755,
 0.847824,
 0.854258,
 0.874916,
 0.918058,
 0.934924,
 0.985091,
 0.994563,
 1.07129,
 1.26215,
 1.26789,
 1.29716,
 1.31747,
 1.36306,
 1.38795,
 1.38952,
 1.40596,
 1.47764,
 1.56872,
 1.58318,
 1.63318,
 1.67675,
 1.84704,
 2.00907,
 2.11734,
 2.24543,
 2.42027,
 2.43474,
 2.56124,
 2.59741,
 2.65974,
 2.66459,
 2.79695,
 2.83653,
 2.89103,
 3.14682,
 3.46321,
 3.60174,
 3.77805,
 4.40907,
 4.53845,
 4.61761,
 5.11489,
 5.62144,
 6.46466,
 7.38818,
 7.75759,
 8.3117,
 9.23522,
 9.23522
])

Dth = np.array([
 0.000569779,
 0.000390302,
 0.000273294,
 0.00057369,
 0.000266974,
 0.000767594,
 0.00039022,
 0.000411818,
 0.000303639,
 0.000597545,
 0.000598725,
 0.000446833,
 0.000390869,
 0.000696567,
 0.000929224,
 0.000307151,
 0.00134012,
 0.000409959,
 0.000864726,
 0.000335311,
 0.00107666,
 0.000595417,
 0.000342933,
 0.00185845,
 0.000385768,
 0.000938037,
 0.000602111,
 0.00068618,
 0.000611549,
 0.000532878,
 0.000484413,
 0.000459396,
 0.00131337,
 0.0009393,
 0.000271925,
 0.000928965,
 0.00080831,
 0.00054163,
 0.000530286,
 0.00128478,
 0.000883429,
 0.0013043,
 0.00261553,
 0.000685197,
 0.00139936,
 0.000876801,
 0.00154184,
 0.0012353,
 0.00195058,
 0.00127527,
 0.00377708,
 0.00197341,
 0.00250896,
 0.00376758,
 0.00607149,
 0.00317883,
 0.00293667,
 0.00328417,
 0.00573614,
 0.00420685,
 0.00546812,
 0.00662834,
 0.006916,
 0.0124427,
 0.00796102,
 0.013219,
 0.0137269,
 0.0150224,
 0.0174393,
 0.027168,
 0.0352153,
 0.0347464,
 0.0393012,
 0.031808
])

"""

#data_run_14 alpha = 0.5
vphi = np.array([
 0.0152451,
 0.0617725,
 0.0932644,
 0.1057,
 0.169572,
 0.177964,
 0.19263,
 0.248705,
 0.261683,
 0.274412,
 0.285389,
 0.385424,
 0.435823,
 0.437378,
 0.45503,
 0.518526,
 0.607615,
 0.628353,
 0.630259,
 0.671109,
 0.75578,
 0.776569,
 0.855332,
 0.928875,
 0.932644,
 0.988602,
 1.00501,
 1.0473,
 1.10733,
 1.14417,
 1.14841,
 1.17327,
 1.17526,
 1.22466,
 1.2256,
 1.23109,
 1.25406,
 1.32802,
 1.40255,
 1.48793,
 1.57664,
 1.5969,
 1.63354,
 1.86679,
 1.99972,
 2.0829,
 2.26105,
 2.41428,
 2.48937,
 2.85737,
 2.97764,
 3.04067,
 3.04903,
 3.26159,
 3.32694,
 3.55369,
 3.55535,
 3.68061,
 3.81693,
 4.01227,
 4.18301,
 4.29845,
 4.45448,
 4.59793,
 4.83993,
 5.07358,
 5.45014,
 6.00223,
 6.94597,
 8.32384,
 9.51296,
 10.1472,
 11.8912
])

Dth = np.array([
 0.000378523,
 0.000600854,
 0.000766539,
 0.000602985,
 0.000399879,
 0.000395364,
 0.000299973,
 0.000347517,
 0.000354077,
 0.00112859,
 0.000406265,
 0.000512452,
 0.000335422,
 0.000470627,
 0.000283534,
 0.000558156,
 0.000469081,
 0.000363464,
 0.000520882,
 0.000382203,
 0.000380108,
 0.000584498,
 0.000614917,
 0.000381472,
 0.000807847,
 0.000442653,
 0.000462552,
 0.00116272,
 0.000969322,
 0.00143537,
 0.000454668,
 0.00143512,
 0.000619769,
 0.000647513,
 0.000513785,
 0.00190257,
 0.00045509,
 0.000934541,
 0.000643639,
 0.000761034,
 0.00048402,
 0.000757781,
 0.000733456,
 0.000834376,
 0.00125227,
 0.0013343,
 0.000875601,
 0.00138782,
 0.00119815,
 0.00164189,
 0.0021593,
 0.00278223,
 0.00217065,
 0.00344058,
 0.0049143,
 0.00422018,
 0.0032302,
 0.00436915,
 0.00379292,
 0.00487055,
 0.00658683,
 0.00606215,
 0.00638879,
 0.007289,
 0.00929521,
 0.0102282,
 0.0130889,
 0.0140999,
 0.0211804,
 0.0336294,
 0.0446268,
 0.0331042,
 0.033603
])

#"""
V0 = 10

##""""
##
####

plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)
#plt.plot(v_centers, D_var, '-o', label=r"$D_{0}(v)$")
plt.plot(v_centers, D_msd, '-o', label=r"$D_{1}(v)$")
plt.plot(v_centers, D_sigma, '-s', label=r"$D_{2}(v)$")
plt.scatter(vphi, Dth,s=50, color='red', marker='o', label=r"$D_{QLT}$")
plt.xlim([v_centers[0], v_centers[-1]])
#plt.ylim([0, 0.006])
#plt.axvline(v_phase, color='k', ls='--', label=r'$v_{ph}$')
plt.xlabel(r"$v$")
plt.ylabel(r"$D(v)$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_fig, "diffusion_coefficient_vs_velocity.pdf"), dpi=dpi)
plt.show()


plt.figure()
for ib in range(0, nv_bins, 4):
    plt.plot(tau_phys, sigma[ib]**2, alpha=0.4)
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\sigma^2(\Delta v)$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

#DIFFUSION EXPONENT α(v):  σ² ~ τ^α

n_v = np.zeros(nv_bins)

for ib in range(nv_bins):
    y = sigma[ib, :Nfit]**2
    n_v[ib] = np.polyfit(np.log(tau_phys[:Nfit]), np.log(y), 1)[0]

plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)
plt.plot(v_centers, n_v, '--',color = "blue")
plt.scatter(v_centers, n_v, color='blue', s=40)
#plt.axhline(1.0, ls='--', color='k', label='Normal diffusion')
#plt.axvline(v_phase, color='r', ls=':')
plt.xlabel(r"$v$")
plt.ylabel(r"$n$  ($\sigma^2 \sim \tau^{\,n}$)")
#plt.title(r"$\sigma^2 \sim \tau^\alpha$")
#plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_fig, "diffusion_exponent_vs_velocity.pdf"), dpi=dpi)
plt.show()


# NON-GAUSSIANITY: KURTOSIS κ(v, τ)
kurt = np.zeros((nv_bins, n_tau))

for ib in range(nv_bins):
    dv = pdf_centers[f"bin_{ib}"]
    for it in range(n_tau):
        P = all_pdfs[ib, it]
        P /= np.trapz(P, dv)
        mean = np.trapz(dv * P, dv)
        var  = np.trapz((dv - mean)**2 * P, dv)
        kurt[ib, it] = np.trapz((dv - mean)**4 * P, dv) / var**2

plt.figure()
plt.imshow(kurt, aspect='auto', origin='lower',interpolation="bilinear",extent=[tau_phys[0], tau_phys[-1],v_centers[0], v_centers[-1]])
plt.colorbar(label="Kurtosis κ")
#plt.axhline(v_phase, color='w', ls='--')
plt.xlabel(r"$\tau$")
plt.ylabel(r"$v$")
plt.title("Non-Gaussianity map")
plt.tight_layout()
plt.show()



# Select velocity bin (e.g. closest to phase velocity)
ib = np.argmin(np.abs(v_centers - V0))
print("bin vel close to vphase",v_centers[ib])
dv = pdf_centers[f"bin_{ib}"]

# Choose tau indices to plot
tau_list = [1, 4, 8, 12, 16, 20]

plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)

for it in tau_list:
    P = all_pdfs[ib, it]

    # Normalize again for safety (optional but recommended)
    norm = np.trapz(P, dv)
    if norm > 0:
        P = P / norm

    plt.plot(dv, P, lw=1.8, label=rf"$\tau = {it:.2f}$")

plt.xlabel(r"$\Delta v$")
plt.ylabel(r"$P(\Delta v,\tau)$")
plt.yscale('log')
#plt.title(rf"PDF of velocity increments at $v_0 \approx {v_centers[ib]:.2f}$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_fig, "PDF_vs_tau.pdf"), dpi=dpi)
plt.show()



plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)
plt.imshow(all_var, aspect='auto', origin='lower',interpolation= "bilinear",extent=[tau_phys[0], tau_phys[-1],v_centers[0],v_centers[-1]],cmap='RdBu_r')#,vmin=-0.5,vmax=0.5)
plt.colorbar(label=r"$\sigma^2$")
plt.xlabel(r"$\tau$")
plt.ylabel(r"$v$")
#plt.title("Velocity diffusion map")
#plt.tight_layout()
plt.savefig(pjoin(path_fig, "sigma_map.pdf"), dpi=dpi)
plt.show()



#test plot for a single velocity bin
#ib = 10
ib = np.argmin(np.abs(v_centers - v_phase))

msd_array = all_var[ib]                 # <(Δv)^2> from MSD
sigma2 = all_sigma_fit[ib]**2           # σ^2 from Gaussian fit

# ---------- Linear fits: sigma^2 ~ a*tau + b ----------
slope_msd, intercept_msd = np.polyfit(tau_phys[:Nfit], msd_array[:Nfit], 1)

slope_sigma, intercept_sigma = np.polyfit(tau_phys[:Nfit], sigma2[:Nfit], 1)

# Fitted straight lines
sigma2_msd_fit   = slope_msd   * tau_phys[:Nfit] + intercept_msd
sigma2_gauss_fit = slope_sigma * tau_phys[:Nfit] + intercept_sigma

# ---------- Plot ----------
plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)

# Scatter points
plt.scatter(tau_phys[:Nfit], msd_array[:Nfit],s=30, color='blue', label=r'$\sigma^2_{1}$')

plt.scatter(tau_phys[:Nfit], sigma2[:Nfit],s=30, color='red', label=r'$\sigma^2_{2}$')

# Fitted lines
plt.plot(tau_phys[:Nfit], sigma2_msd_fit,color='blue', linewidth=2)
#plt.plot(tau_steps[:Nfit], sigma2_msd_fit,color='blue', linewidth=2)

plt.plot(tau_phys[:Nfit], sigma2_gauss_fit,color='red', linewidth=2)

plt.xlabel(r'$\tau$')
plt.ylabel(r'$\sigma^2$')
#plt.title(rf"Diffusion fit at $v \approx {v_centers[ib]:.2f}$")
plt.grid(alpha=0.35)
plt.legend()
plt.tight_layout()

plt.savefig(pjoin(path_fig, f"sigma2_vs_tau_bin{ib}.pdf"), dpi=900)
plt.show()

###$$$


# ==========================================================
# SELF-SIMILARITY TEST: PDF COLLAPSE
# ==========================================================

# Use velocity bin closest to phase velocity
ib = np.argmin(np.abs(v_centers - V0))
dv = pdf_centers[f"bin_{ib}"]

# Choose τ values (spread across range)
tau_list = [1, 4, 8, 12, 16, 20]

plt.figure(figsize= figsize / 25.4, constrained_layout=True, dpi=ppi)

for it in tau_list:
    P = all_pdfs[ib, it]
    s = sigma[ib, it]

    #norm = np.trapz(P, dv)
    #if norm > 0:
    #    P = P / norm
    #    print("norm:", norm)

    # Self-similar rescaling
    #alpha = 0.9
    #scale = tau_phys[it]**alpha
    #plt.plot(dv / scale, P * scale, lw=1.8, alpha=0.7, label=rf"$\tau = {tau_phys[it]:.2f}$")
    plt.plot(dv / s, P * s, lw=1.8, alpha=0.7, label=rf"$\tau = {it:.2f}$")

# Reference Gaussian
x = np.linspace(-5, 5, 400)
plt.plot(x, np.exp(-x**2/2)/np.sqrt(2*np.pi),'k--', lw=2, label="Gaussian")

plt.xlabel(r"$\Delta v / \sigma$")
plt.ylabel(r"$\sigma\, P(\Delta v,\tau)$")
plt.yscale('log')
plt.xlim(-10,10)
#plt.title(rf"PDF self-similarity test at $v_0 \approx {v_centers[ib]:.2f}$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(pjoin(path_fig, "PDF_self_similarity_test.pdf"), dpi=dpi)
plt.show()
