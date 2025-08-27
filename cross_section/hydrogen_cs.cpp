#include "hydrogen_cs.h"
// Maximum number of data points in .dat file
#define MAX_POINTS 500

// Arrays to store cross-section data
double elastic_energies[MAX_POINTS];
double elastic_cross_sections[MAX_POINTS];
double excitation_energies[MAX_POINTS];
double excitation_cross_sections[MAX_POINTS];
double ionization_energies[MAX_POINTS];
double ionization_cross_sections[MAX_POINTS];
int elastic_count = 0;
int excitation_count = 0;
int ionization_count = 0;

// Load data from a .dat file
void load_file(const char* filename, double energies[], double cross_sections[], int* count)
{
    FILE* file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open file !!! %s\n", filename);
        *count = 0;
        return;
    }

    double energy, cross_section;
    *count = 0;
    char line[100];

    while (fgets(line, sizeof(line), file) && *count < MAX_POINTS)
    {
        if (line[0] == '#') continue;
        if (sscanf(line, "%lf %lf", &energy, &cross_section) == 2)
        {
            energies[*count] = energy;
            cross_sections[*count] = cross_section;
            (*count)++;
        }
    }

    fclose(file);
    if (*count == 0)
    {
        printf("Error: No data in file %s\n", filename);
    }
}

// Interpolate cross-section value for a given energy
double interpolate(double energy, double energies[], double cross_sections[], int count)
{
    if (count == 0 || energy < energies[0] || energy > energies[count - 1])
    {
        return 0.0;
    }

    for (int i = 0; i < count - 1; i++)
    {
        if (energy >= energies[i] && energy <= energies[i + 1])
        {
            double e1 = energies[i];
            double e2 = energies[i + 1];
            double cs1 = cross_sections[i];
            double cs2 = cross_sections[i + 1];
            return cs1 + (cs2 - cs1) * (energy - e1) / (e2 - e1);
        }
    }

    return 0.0;
}

// Set up Hydrogen cross-section form .dat files
void setup_hydrogen_cross_sections(const char* data_path)
{
    char filename[100];
    snprintf(filename, sizeof(filename), "%s/hydrogen_elastic_cs.dat", data_path);
    load_file(filename, elastic_energies, elastic_cross_sections, &elastic_count);

    snprintf(filename, sizeof(filename), "%s/hydrogen_excitation_cs.dat", data_path);
    load_file(filename, excitation_energies, excitation_cross_sections, &excitation_count);

    snprintf(filename, sizeof(filename), "%s/hydrogen_ionization_cs.dat", data_path);
    load_file(filename, ionization_energies, ionization_cross_sections, &ionization_count);
}

// Hydrogen cross-section functions
double compute_elastic_CS_h(double energy, Domain &domain)
{
    if (!domain.enable_elastic_collision) return 0.0;
    return interpolate(energy, elastic_energies, elastic_cross_sections, elastic_count);
}

double compute_excitation_CS_h(double energy, Domain &domain)
{
    if (!domain.enable_excitation_collision || energy < E_EXC_TH_H) return 0.0;
    return interpolate(energy, excitation_energies, excitation_cross_sections, excitation_count);
}

double compute_ionization_CS_h(double energy, Domain &domain)
{
    if (!domain.enable_ionization_collision || energy < E_ION_TH_H) return 0.0;
    return interpolate(energy, ionization_energies, ionization_cross_sections, ionization_count);
}


double compute_mex_cs_pion(double energy, Domain &domain)
{
    if (!domain.enable_pion_elastic) return 0.0;

    return (1e-20 * 5.74 * pow((energy / 13.61), -0.5765)) / (1 + pow((energy / 27.9), 1.1605));

}

double compute_e_detach(double energy, Domain &domain)
{
    if (!domain.enable_e_detach_collision || energy <= 2.25) return 0.0;
    return ((1e-20 * 4.19e-2 * pow((energy - 2.25) / 13.61, 1.89)) / (1.0 + pow((energy - 2.25) / 178.0, 1.66) 
    + pow((energy - 2.25) / 1040.0, 2.76)) + (1e-20 * 16.5 * pow((energy - 2.25) / 13.61, 1.088)) / (1.0 + pow((energy - 2.25) / 5.33, 1.254)));
    
}