#include "init.h"

Init::Init(Species &species, Domain &domain) : species(species),domain(domain)
{

    auto initilization1 = INIParser::loadtypeextract(species.initialization_pos);
    auto initilization2 = INIParser::loadtypeextract(species.initialization_vel);

    auto [init_type1, n1, amplitude1] = initilization1;
    auto [init_type2, n2, amplitude2] = initilization2;

    //display::print(init_type1);
    //display::print(init_type2, n2, amplitude2);
    double k_launch = 0.0;

    for (int p = 0; p < species.numparticle; p++) 
    {
        double x = 0;
        double vx = 0, vy = 0, vz = 0;

        // ----- Position initialization -----
        if (init_type1 == "random") 
        {
            x = domain.x0 + domain.xL * rnd();
        }
        
        else if (init_type1 == "uniform") 
        {
            x = domain.x0 + p * (domain.xL / (species.numparticle - 1));
        }
        else if (init_type1 == "sin" || init_type1 == "cos") 
        {
            double x0 = domain.x0 + p * (domain.xL / (species.numparticle - 1));
            double k = 2 * Const::PI * n1 / domain.xL;
            double A = amplitude1;

            if (init_type1 == "sin")
            {
                x = x0 + A * sin(k * x0);
            } 
            else
            {
                x = x0 + A * cos(k * x0);
            }
        }
        else if (init_type1 == "extend")
        {
            int start = static_cast<int>(n1);
            int end = static_cast<int>(amplitude1);
            x = start + (end - start) * rnd();  // uniform random in [start, end]
        }

        //out of bound case
        if (x >= domain.xL)
        {
            x -= domain.xL;
        }
        else if (x < 0)
        {
            x += domain.xL;
        }
        // ----- Velocity initialization -----

        vx = SampleVel(species) + species.vs * domain.vel_norm;

        if (init_type2 == "sin" || init_type2 == "cos") 
        {
            double k = 2 * Const::PI * n2 / domain.xL;
            k_launch = k;
            
            double A = amplitude2 * domain.vel_norm;

            if (init_type2 == "sin")
            {
                vx += A * sin(k * x);  // NOTE: use x, not x0
            }   
            else
            {
                vx += A * cos(k * x);
            }    
        }
        else
        {
            display::print("ERROR! Wrong velcoity init");
            exit(-1);
        }
        // Normalize velocity
        vx /= domain.vel_norm;
        vy = 0.0;
        vz = 0.0;

        species.AddParticle(Particle(x, vx, vy, vz, 0));
    }
    display::print("launched wave k:", k_launch);
}

double Init::SampleVel(Species &species)
{
    //double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);
    double v_th =  sqrt(2*Const::K_b*species.temp*Const::EV_to_K/species.defaultmass);
    double vt = v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5) ;//+ domain.v_i*domain.wp * domain.LD;
    //double vt = v_th * sqrt(2) * (rnd()*rnd()*rnd() - 1.5) ;
    return vt;
}

double Init::SampleVel(Species &species, double temp)
{
    //double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);
    double v_th =  sqrt(2*Const::K_b*temp*Const::EV_to_K/species.defaultmass);
    double vt = v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5); //+ domain.v_i*domain.wp * domain.LD;
    return vt;
}
