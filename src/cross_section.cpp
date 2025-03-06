#include "cross_sections.h"

// Elastic cross-section  for e- / Ar collision
double compute_elastic_CS(double energy, Domain &domain)
{ 
    if (!domain.enable_elastic_collision) return 0.0;

    return (1e-20 * fabs(6.0 / pow(1.0 + (energy / 0.1) + pow(energy / 0.6, 2.0), 3.3) -
    1.1 * pow(energy, 1.4) / (1.0 + pow(energy / 15.0, 1.2)) / sqrt(1.0 + pow(energy / 5.5, 2.5)
     + pow(energy / 60.0, 4.1)) + 0.05 / pow(1.0 + energy / 10.0, 2.0) + 0.01 * pow(energy, 3.0) / (1.0 + pow(energy / 12.0, 6.0))));
}

// Excitation cross-section  for e- / Ar collision
double compute_excitation_CS(double energy, Domain &domain)
{
    
    if (!domain.enable_excitation_collision) return 0.0;

    if (energy > E_EXC_TH)
    {
        return (1e-20 * (0.034 * pow(energy - E_EXC_TH, 1.1) * (1.0 + pow(energy / 15.0, 2.8)) / (1.0 + pow(energy / 23.0, 5.5))
            + 0.023 * (energy - E_EXC_TH) / pow(1.0 + energy / 80.0, 1.9)));
    }
    return 0.0;
}

// Ionization cross-section for e- / Ar collision
double compute_ionization_CS(double energy, Domain &domain)
{
    
    if (!domain.enable_ionization_collision) return 0.0;

    if (energy > E_ION_TH)
    {
        return (1e-20 * (970.0 * (energy - E_ION_TH) / pow(70.0 + energy, 2.0)
            + 0.06 * pow(energy - E_ION_TH, 2.0) * exp(-energy / 9)));
    }
    return 0.0;
}
