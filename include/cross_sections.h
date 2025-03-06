#ifndef _CROSS_SECTIONS_H_
#define _CROSS_SECTIONS_H_

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
double compute_elastic_CS(double energy, Domain &domain);      // Elastic cross-section
double compute_excitation_CS(double energy, Domain &domain);   // Excitation cross-section
double compute_ionization_CS(double energy, Domain &domain);   // Ionization cross-section

#endif  // _CROSS_SECTIONS_H_
