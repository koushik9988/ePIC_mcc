#ifndef _SPECIES_H_
#define _SPECIES_H_

#include <iostream>
#include <list>
#include <vector>

#include <thread>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <mutex>
#include <functional>
#include <init.h>
#include <tuple>
#include "slap.h"
#include "domain.h"
#include "slap.h"

using namespace std;

class Domain;


/**
 * @class Particle
 * @brief Particle class to hold particel position and velocity data.
 */
/// @brief
class Particle
{
    public:
    double x;
    //double vel[3];
    double vx, vy, vz;
    double mass;
    double spwt; //specific weight
    //store velocity at (t-0.5dt)
    double pvx,pvy,pvz;

    int id; //particle id
    //particle constructor
    Particle(double x, double vx, double vy, double vz, int id, double mass = NAN , double spwt = NAN)
    {
        this->x = x;
        this->vx = vx;
        this->vy = vy;
        this->vz = vz;
        this->pvx = 0;
        this->pvy = 0;
        this->pvz = 0;
        this->id = id; //default id
        this->mass = mass;
        this->spwt = spwt;
    }

};

/**
 * @class Species
 * @brief represent indivitual species and each species propery and various methods .
 */

/// @brief
class Species
{
    public:
    /// species name
    string name;
    /// species mass
    double defaultmass;
    /// species charge
    double charge;
    /// @brief specific weight
    double defaultspwt;
    /// @brief species temparature
    double temp;
    /// @brief density
    vec<double> den;
    vec<double> velmesh;
    /// @brief no of simulation particle
    int numparticle;
    /// @brief 
    int charge_sig;

    //collison rate
    vec<double>coll_rate;
    
    std::string initialization_pos;
    std::string initialization_vel;

    double vs;
    double fract_den;

    // the "part_list" is a linked list that holds instances of the Particle class.
    //list/vector is a template class and Particle is template argument(similar to regular int,double etc)
    //specifying type of element it will hold which is instaances of "Particle" class here. 
    //(C++ STL containers such as vector,list can hold generic data type like class etc.)
    //list<Particle> part_list;
    //or
    
    vector<Particle> part_list;

    vec<double> velprev;


    //vector<vec<double>> buffers;
    //buffers.reserve();

    /**
     * @brief Constructor for species class.
     * @param mass mass of species.
     * @param charge charge of species
     * @param spwt specific weight of species
     * @param temp temparature of species
     * @param numparticle number pf simulation particle
     */
    Species(string name, double defaultmass, double charge, double defaultspwt, double temp, int numparticle, double vs, double fract_den, std:: string initialization_pos, std:: string initialization_vel, Domain &domain);

    //declare member functions or methods
    /**
     * @brief function to Add Particle.
     * @param part particle class instances.
     */
    void AddParticle(Particle part);
    /**
     * @brief function to move the particle.
     * @param species_list list of species.
     * @param sub_cycle sub cycle interval.
     */
    void Push_species(vector<Species> &species_list, int sub_cycle);
    /**
     * @brief helper function for particle mover parallization.
     * @param start staring index for one thread.
     * @param end end index for one thread
     * @param sub_cycle sub cycle interval.
     * @param species_list list of species
     */
    void update(int start, int end, int sub_cycle,vector<Species> &species_list);
    /**
     * @brief serial code for particle mover function.
     * @param sub_cycle sub cycle interval.
     * @param species_list list of species
     */
    void Push_species_serial(vector<Species> &species_list, int sub_cycle);
    /**
     * @brief function to deposit charge on grid points.
     */
    void ScatterSpecies();
    /**
     * @brief function to deposit charge on grid points.
     */
    void ScatterSpecies_serial();
    /**
     * @brief helper function for parallal charge deposition.
     * @param threadidx current thread index 
     * @param start staring index for one thread.
     * @param end end index for one thread
     */
    void parall_deposit(int threadidx, int start, int end);
    /**
     * @brief function of rewind velocity by half time step (which is required for leap-froging).
     */
    void Rewind_species();
    /**
     * @brief function to compute kinetic energy
     * @param species normalizing species 
     */
    //double Compute_KE(Species &species);
    vec<double> Compute_KE(Species &species);
    /**
     * @brief function to compute momentum energy
     * @param species normalizing species 
     */
    //double Compute_Momentum(Species &species);
    std::tuple<double, double, double> Compute_Momentum(Species &species);
    /**
     * @brief function to check if a species is Ion.
     */
    bool IsIon();

    void ScatterVel_serial();

    private:
    Domain &domain;  
};


/**
 * @struct SpeciesParams
 * @brief Structure to hold parameters for initializing a species.
 */
/// @brief
struct SpeciesParams
{
    string name;
    double mass;
    int num;
    double temp;
    int charge_sign;
    double normden;
    double vs;
    string loadtype_pos;
    string loadtype_vel;
};



/**
*@brief function for secondary elctron emission.
* @param species instances of sepcies class.
* @param domain instances of domain class.
*/
void SEE(Species &species, Domain &domain);

#endif 