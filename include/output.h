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
#include "matplotlibcpp.h"

#include "H5Cpp.h"

using namespace H5;
namespace plt = matplotlibcpp;

class Domain;
class Species;

class Output 
{
    public:

    std::map<std::string, Group> particle_groups;
    std::map<std::string, Group> den_subgroups;
    std::map<std::string, Group> vel_subgroups;
    std::map<std::string, Group> collrate_subgroups;

    Output(const std::filesystem::path& outputfolder, Domain& domain);
    //~Output();
    void write_particle_data(int ts, Species& species);
    //void write_particle_data(H5::Group& group, int ts, Species& species);
    void write_den_data(int ts,  Species& species);
    void write_collrate_data(int ts, Species& species);
    void write_vel_data(int ts,  Species& species);
    void write_field_data(int ts);
    //void write_ke(int ts,std::vector<Species> &species_list);
    void write_ke();
    void write_m();
    void storeKE_to_matrix(int ts, std::vector<Species> &species_list);
    void storem_to_matrix(int ts, std::vector<Species> &species_list);
    void printmatrix(int row, int col, double **matrix);
    void diagnostics(int ts, std::vector<Species> &species_list);
    void write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT,double density, int save_fig, int normscheme, 
        int subcycleint,double LDe, double LDi, double wpe, double wpi,int spno, double GAS_DENSITY, double max_electron_collision_freq);
    void write_species_metadata(std::vector<Species>& species_list);

    void write_avg_collision_freq(int ts);
    void write_alpha_vs_time(int ts);

    void write_extra(double z1, double z2, double z3);

    //data structure to temporarily store kinetic energy and momentum data.

    Matrix<double> store_ke;
    Matrix<double> store_m;

    int sp_no ;//= species_list.size();
    int t_step;// = int(domain.NUM_TS/domain.write_interval) + 1 ;
    
    int precision;

    //plotting flag
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
    H5File file; // Declare H5::H5File object to handle HDF5 file operations
    Group field_data_group;
    Group time_group;
    Group metadata_group;
    Group metadata_species;
    //std::ofstream report;
    //std::vector<Species> species_list;
   
};

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
