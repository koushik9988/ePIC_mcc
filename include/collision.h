#ifndef _COLLISIONAL_MODEL_H_
#define _COLLISIONAL_MODEL_H_

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "domain.h"
#include "species.h"
#include "init.h"
#include <tuple>
#include "slap.h"

#include "Ar_cs.h"
#include "hydrogen_cs.h"


// Define gas types
#define ARGON 0
#define HYDROGEN 1

constexpr int CS_RANGES = 2000000; //CS_RANGES is a compile-time constant
constexpr double DE_CS = 0.005;  
//using cross_section = std::array<double, CS_RANGES>;
using cross_section = std::vector<double>;

class CollisionHandler
{
    private:
    Domain &domain;

    //total cross section
    cross_section sigma_tot_e;
    cross_section sigma_tot_pion;
    cross_section sigma_tot_nion;

    //indivisual cross ection
    cross_section sigma_ela;
    cross_section sigma_exc;
    cross_section sigma_ionz;
    cross_section sigma_mex_pion;
    cross_section sigma_e_detach;
    
    public:
    int gas;
    CollisionHandler(Domain& domain, int gas, const char* data_path);
    CollisionHandler(Domain &domain, int gas);
    void set_electron_cross_sections();
    void set_pion_cross_sections();
    void calc_total_cross_sections();
    void collision_electron(double xe, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2);
    void collision_ion(double xe, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2, int eindex, Species &species1, Species &species2);
    void collision_detachment(double xe, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2, int eindex, Species &species1, Species &species2, Species &electron_species);
    //void handle_collisions(Species &electron, Species &target_gas);
    void handle_collisions(Species &projectile, Species &target_gas, Species &electron_species);
    double max_electron_coll_freq();
    //double average_collision_frequency(Species &electron);
    double average_collision_frequency(Species &projectile);
    std::pair<double, double> electron_min_max_mean_free_path();
    void coll_rate(Species &projectile, Species &target_gas);
};

#endif  // _COLLISIONAL_MODEL_H_
