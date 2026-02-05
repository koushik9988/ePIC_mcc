#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <fstream>
#include <vector>
#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <map>
#include <memory>
#include "domain.h"
#include "species.h"
#include "math.h"
#include "slap.h"
#include "H5Cpp.h"


#ifdef ENABLE_PLOTTING
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

using namespace H5;

class Domain;
class Species;
class CollisionHandler;


/**
 * @class Output
 * @brief Manages output operations for the simulation, including writing particle data, field data,
 *        and diagnostics to HDF5 files, as well as generating plots if enabled.    
 */


/// @brief  
class Output 
{
    public:
    /// 
    std::map<std::string, Group> particle_groups;
    std::map<std::string, Group> den_subgroups;
    std::map<std::string, Group> vel_subgroups;
    std::map<std::string, Group> collrate_subgroups;

    /**
     * @brief Constructor for the Output class.
     * @param outputfolder Path to the output directory.
     * @param domain Reference to the simulation domain.
     * @note Initializes HDF5 file and groups for storing simulation data.
     */
    Output(const std::filesystem::path& outputfolder, Domain& domain);
    //~Output();
    /**
     * @brief function to write particle data (phase space data) to HDF5 file.
     * @param ts Current time step.
     * @param species Reference to the species whose data is to be written.
     */
    void write_particle_data(int ts, Species& species);
    /**
     * @brief function to write density data to HDF5 file.
     * @param ts Current time step.
     * @param species Reference to the species whose density data is to be written.
     */
    void write_den_data(int ts,  Species& species);
    /**
     * @brief function to write collision rate.
     * @param ts Current time step.
     * @param species Reference to the species whose velocity data is to be written.
     */
    void write_collrate_data(int ts, Species& species);
    /**
     * @brief function to write velocity data to HDF5 file.
     * @param ts Current time step.
     * @param species Reference to the species whose velocity data is to be written.
     */
    void write_vel_data(int ts,  Species& species);
    /**
     * @brief function to write field data (potential, electric field, charge density) to HDF5 file.
     * @param ts Current time step.
     */
    void write_field_data(int ts);
    /**
     * @brief function to write kinetic energy and momentum data to HDF5 file.
     */
    void write_ke();
    /**
     * @brief function to write momentum data to HDF5 file.
     */
    void write_m();
    /**
     * @brief function to store kinetic energy data into a matrix for later writing.
     * @param ts Current time step.
     * @param species_list List of species in the simulation.
     */
    void storeKE_to_matrix(int ts, std::vector<Species> &species_list);
    /**
     * @brief function to store momentum data into a matrix for later writing.
     * @param ts Current time step.
     * @param species_list List of species in the simulation.
     */
    void storem_to_matrix(int ts, std::vector<Species> &species_list);
    void printmatrix(int row, int col, double **matrix);
    /**
     * @brief function to write diagnostics data to putput screen
     * @param ts Current time step.
     * @param species_list List of species in the simulation.
     */
    void diagnostics(int ts, std::vector<Species> &species_list);
    /**
    *@brief function to write simulation metadata to HDF5 file.
     * @param NC number of grid points.
     * @param NUM_TS Total simulation time step.
     * @param write_int data writing interval.
     * @param write_int_phase phase space data writing interval.
     * @param DT time step.
     * @param density plasma density.
     * @param save_fig flag to save plots.
     * @param normscheme normalization scheme flag.
     * @param subcycleint ion sub-cycling interval.
     * @param LDe electron debye length.
     * @param LDi ion debye length.
     * @param wpe electron plasma frequency.
     * @param wpi ion plasma frequency.
     * @param spno number of species.
     * @param GAS_DENSITY neutral gas density.
     * @param max_electron_collision_freq maximum electron collision frequency.
     */
    void write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT,double density, int save_fig, int normscheme, 
        int subcycleint,double LDe, double LDi, double wpe, double wpi,int spno, double GAS_DENSITY, double max_electron_collision_freq, double B, double theta, double azimuth);
    /**
     * @brief function to write species metadata to HDF5 file.
     * @param species_list List of species in the simulation.
     */
    void write_species_metadata(std::vector<Species>& species_list);

    void write_avg_collision_freq(int ts);
    void write_alpha_vs_time(int ts);

    void write_extra(double z1, double z2, double z3);

    void store_avg_coll_freq_to_matrix(int ts, Species &species, double avg_coll_freq, std::vector<Species> &species_list);
    void write_avg_coll_freq();

    void write_dump_file(const std::filesystem::path& dump_path,int ts,std::vector<Species>& species_list);


    ///data structure to temporarily store kinetic energy and momentum data.
    Matrix<double> store_ke;
    Matrix<double> store_m;
    Matrix<double> store_avg_coll_freq;

    int sp_no ;//= species_list.size();
    int t_step;// = int(domain.NUM_TS/domain.write_interval) + 1 ;
    
    int precision;

    /// @brief flags to control what data to plot and save
    int Energy_plot;
    int Potentialfield_plot;
    int Chargedensity_plot;
    int keflag;
    int peflag;
    int teflag;
    int phase_plot;
    int dft_flag;
    int species_index;
    int ke_components;
    int coll_freq_plot;

    /// @brief data structures to hold time series data for diagnostics    
    std::vector<double> time_steps;
    std::vector<double> kinetic_energy;
    std::vector<double> potential_energy;
    std::vector<double> total_energy;
    std::vector<double>Ke_x;
    std::vector<double>Ke_y;
    std::vector<double>Ke_z;
    std::vector<double> nu_avg;

    private:
    //Species &species;
    std::filesystem::path outputfolder;
    Domain& domain;
    CollisionHandler *coll;
    H5File file; // Declare H5::H5File object to handle HDF5 file operations
    Group field_data_group;
    Group time_group;
    Group metadata_group;
    Group metadata_species;
    //std::ofstream report;
    //std::vector<Species> species_list;
   
};


/**
 * @namespace display
 * @brief Utility namespace for displaying messages to the console.
 */
namespace display
{
    template<typename T>
    void print(const T& value) 
    {
        std::cout << value << std::endl;
    }
    // Recursive template function to print multiple arguments
    template<typename T, typename... Args>
    void print(const T& value, Args&&... args) 
    {
        std::cout << value;
        print(std::forward<Args>(args)...); // Recursive call with the remaining arguments
    }

}

#endif
