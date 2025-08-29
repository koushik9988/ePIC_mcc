#ifndef _INIT_H_
#define _INIT_H_

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "species.h"
#include "domain.h"
#include "iniparser.h"
#include <tuple>

class Domain;
class Species;

//using namespace Const;

class Init 
{
    public:
    //in the constructor we pass Species and Domain class instances to access the species and domain class
    //memebr variable and methods 
    Init(Species &species, Domain &domain);   
    static double SampleVel(Species &species, double temp);

    private:
    Species &species;
    Domain &domain;
    double SampleVel(Species &species);
};

#endif 
