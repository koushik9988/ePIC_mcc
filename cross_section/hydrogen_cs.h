#ifndef _HYDROGEN_CROSS_SECTIONS_H_
#define _HYDROGEN_CROSS_SECTIONS_H_

#include "domain.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Define constants for Hydrogen thresholds
#define E_EXC_TH_H 10.2  // Excitation threshold in eV 
#define E_ION_TH_H 13.6  // Ionization threshold in eV 

// Function to set up Hydrogen cross-section data
void setup_hydrogen_cross_sections(const char* data_path);

//Functions to calculate Hydrogen cross-sections
double compute_elastic_CS_h(double energy, Domain &domain);
double compute_excitation_CS_h(double energy, Domain &domain);
double compute_ionization_CS_h(double energy, Domain &domain);

//H^+ + H2 -> elastic cross section
double compute_mex_cs_pion(double energy, Domain &domain);
double compute_mex_cs_nion(double energy, Domain &domain);
double compute_e_detach(double energy, Domain &domain);

#endif  // _HYDROGEN_CROSS_SECTIONS_H_