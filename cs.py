import numpy as np
import matplotlib.pyplot as plt

def sigma_Hp_H2_momentum_transfer(E_eV):
    sigma0 = 1e-16  # cm²
    a1 = 5.74
    a2 = -0.5765
    a3 = 27.9   # eV
    a4 = 1.737
    ER = 13.61  # eV

    numerator = sigma0 * a1 * (E_eV / ER)**a2
    denominator = 1 + (E_eV / a3)**(a2 + a4)
    return numerator / denominator

def sigma_hminus_h2_momentum_transfer(E_eV):
    sigma0 = 1e-16  # cm²
    ER = 13.61      # eV
    a1 = 6.36
    a2 = -0.337
    a3 = 5.5        # eV
    a4 = 0.83
    a5 = 26.0       # eV
    a6 = 1.766

    numerator = sigma0 * a1 * (E_eV / ER)**a2
    denominator = 1 + (E_eV / a3)**(a2 + a4) + (E_eV / a5)**(a2 + a6)
    return numerator / denominator

E_vals = np.logspace(-2, 6, 10000)

# Calculate cross sections
sigma_Hp = sigma_Hp_H2_momentum_transfer(E_vals)
sigma_Hm = sigma_hminus_h2_momentum_transfer(E_vals)

# Plot
fig, ax = plt.subplots(figsize=(9, 6))
ax.loglog(E_vals, sigma_Hp, label=r'H$^+$ + H$_2$', color='red')
ax.loglog(E_vals, sigma_Hm, label=r'H$^-$ + H$_2$', color='blue')

ax.set_xlabel('Energy (eV)')
ax.set_ylabel(r'Cross Section $\sigma$ (cm$^2$)')
ax.set_title(r'Momentum Transfer / Elastic Cross Sections: H$^\pm$ + H$_2$')
ax.grid(True, which='both', ls='--', alpha=0.5)
ax.legend()
plt.tight_layout()
plt.show()
