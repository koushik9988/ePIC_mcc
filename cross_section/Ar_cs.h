#ifndef _ARGON_CROSS_SECTIONS_H_
#define _ARGON_CROSS_SECTIONS_H_

#include <cmath>
#include <cstdlib>
#include <locale>
#include <ios>         // For std::ios_base
#include <iostream>
#include "domain.h"

class Domain;

// Define constants for excitation and ionization thresholds
constexpr double E_EXC_TH = 11.5;  // Excitation energy threshold in eV
constexpr double E_ION_TH = 15.8;  // Ionization energy threshold in eV

// Function declarations (prototypes) for cross sections
double compute_elastic_CS_ar(double energy, Domain &domain);      // Elastic cross-section
double compute_excitation_CS_ar(double energy, Domain &domain);   // Excitation cross-section
double compute_ionization_CS_ar(double energy, Domain &domain);   // Ionization cross-section


double compute_iso_CS_ar_ion(double energy, Domain &domain);      // Isotropic elastic scattering cross-section
double compute_back_CS_ar_ion(double energy, Domain &domain);     // Backward elastic scattering cross-section


#endif  // _CROSS_SECTIONS_H_
