#include "init.h"

Init::Init(Species &species, Domain &domain) : species(species),domain(domain)
{

    auto initilization1 = INIParser::loadtypeextract(species.initialization_pos);
    auto initilization2 = INIParser::loadtypeextract(species.initialization_vel);

    auto [init_type1, n1, amplitude1] = initilization1;
    auto [init_type2, n2, amplitude2] = initilization2;

    //display::print(init_type1);
    //display::print(init_type2, n2, amplitude2);
    //double k_launch = 0.0;


    ////test code
    // -----------------------------------------
    bool dump_mode = (init_type1 == "dump" || init_type2 == "dump" || init_type1 == "dump+sin" || init_type1 == "dump+cos");

    std::vector<double> dump_phase;  // x vx vy vz id
    hsize_t dump_np = 0;

    if (dump_mode)
    {
        H5File file("dump.h5", H5F_ACC_RDONLY);

        std::string dset_path = "/dump/particles/" + species.name + "/phase";
        DataSet dset = file.openDataSet(dset_path);

        DataSpace space = dset.getSpace();
        hsize_t dims[2];
        space.getSimpleExtentDims(dims);

        dump_np = dims[0];

        if (dump_np != static_cast<hsize_t>(species.numparticle))
        {
            //std::cerr << "ERROR: dump particle count mismatch for species "<< species.name << std::endl;
            //exit(-1);
        }

        dump_phase.resize(dump_np * 5);
        dset.read(dump_phase.data(), PredType::NATIVE_DOUBLE);

        file.close();

        //display::print("Restarting sim species:", species.name," with particles:", dump_np);
        // Overwrite particle count from dump
        species.numparticle = dump_np;
        display::print("Loaded particles from dump file.",dump_np);

    }

    //////

    for (int p = 0; p < species.numparticle; p++) 
    {
        double x = 0;
        double vx = 0, vy = 0, vz = 0;

        int idx = 5 * p;

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
            
            //acceptance-rejection sampling (may better reproduce density)
            double k = 2 * Const::PI * n1 / domain.xL;
            double A = amplitude1;

            double pdf_env = 1/(1 + A);

            while(true)
            {
                double x_try = domain.x0 + domain.xL * rnd();

                double s = (init_type1 == "sin") ? sin(k * x_try) : cos(k * x_try);

                double pdf = (1.0 + A * s) * pdf_env;

                if (rnd() < pdf)
                {
                    x = x_try;
                    break;
                }
            }
            
            
            //double x0 = domain.x0 + p * (domain.xL / (species.numparticle - 1));
            //double k = 2 * Const::PI * n1 / domain.xL;
            //double A = amplitude1;

            //if (init_type1 == "sin")
            //{
            //    x = x0 + A * sin(k * x0);
            //} 
            //else
            //{
            //    x = x0 + A * cos(k * x0);
            //}
        }
        else if (init_type1 == "extend")
        {
            int start = static_cast<int>(n1);
            int end = static_cast<int>(amplitude1);
            x = start + (end - start) * rnd();  // uniform random in [start, end]
        }



        //test code
        ///
        else if (init_type1 == "dump")
        {
            x  = dump_phase[idx + 0];
            //display::print("Loaded particle position:", x);
        }
        else if (init_type1 == "dump+sin" || init_type1 == "dump+cos")
        {
            x  = dump_phase[idx + 0];

            //double k = 2 * Const::PI * n1 / domain.xL;
            double k = 2 * Const::PI * n1 / (domain.xL);
            //k_launch = k;

            double A = amplitude1;

            if (init_type1 == "dump+sin")
            {
                x += A * sin(k * x);
            }   
            else
            {
                x += A * cos(k * x);
            }    
        }

        else if (init_type1 == "tanh+sin")
        {
            
            double x_random;
            //construct density profile
            double RB = domain.x0 + 0.25 * domain.xL;  // left boundary/edge
            double LB = domain.x0 + 0.75 * domain.xL;  // right boundary/edge
            double d  = 0.05 * domain.xL;              // interface width

            while (true)
            {
                // pick a random x and check whether it fit in the profile we need or not /comapre againt a random number 
                x_random = domain.x0 + rnd()*domain.xL;
                double rho = 0.5 * (tanh((x_random - RB)/d)-tanh((x_random - LB)/d)); //multiply by half such that value varies between 0 to 1
                if (rnd() < rho)
                {
                    break;
                }
            }
            ///*
            double A = amplitude1;                 // perturbation amplitude
            //double k = (2.0 * M_PI * 1) / (LB - RB);
            double k =  2 * Const::PI * n1 / (domain.xL);   

            x = x_random; //uncomment#when removed test

            //@if (x > RB && x < LB)
            //@{
                //@x += A * std::sin(k * (x - RB));
                x += A * sin(k*x); //uncomment# when removed test
            //@}

        
            
            //*/
            //x = x_random;
        }

        else if (init_type1 == "tanh")
        {
            
            double x_random;
            //construct density profile
            double RB = domain.x0 + 0.25 * domain.xL;  // left boundary/edge
            double LB = domain.x0 + 0.75 * domain.xL;  // right boundary/edge
            double d  = 0.05 * domain.xL;              // interface width

            while (true)
            {
                // pick a random x and check whether it fit in the profile we need or not /comapre againt a random number 
                x_random = domain.x0 + rnd()*domain.xL;
                double rho = 0.5 * (tanh((x_random - RB)/d)-tanh((x_random - LB)/d)); //multiply by half such that value varies between 0 to 1
                if (rnd() < rho)
                {
                    break;
                }
            }
            x = x_random;
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
            //k_launch = k;
            
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

        else if (init_type2 == "dump")
        {
            vx = dump_phase[idx + 1];
            vy = dump_phase[idx + 2];
            vz = dump_phase[idx + 3];
            //display::print("Loaded particle velocity:", vx, vy, vz);
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
    //display::print("launched wave k:", k_launch);
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
