//implementation file(.cpp file)contains implementation of function and class methods
#include "domain.h"

//instances of class Rnd.
Rnd rnd;

//domain constructor
Domain::Domain(double x0, double dx, int ni):x0(x0), dx(dx), ni(ni)
{
    phi = vec<double>(ni);
    rho = vec<double>(ni);
    ef  = vec<double>(ni);

    collision_rate = vec<double>(ni);

    dft_k = vec<double>((ni/2) +1);
    dft_value = vec<double>((ni/2) +1);
    xL = (ni-1)*dx;
}


void Domain::display(vector<Species> &species)
{
    display::print("\n*********************************************************************************************");
    display::print("         1D Electrostatic Particle-In-Cell Code (ePIC++) with MCC Collision Model ");
    display::print("*********************************************************************************************\n");
    
    display::print(" Simulation Parameters");
    display::print("--------------------------------------------");
    display::print("  Cell Spacing:               ", dx);
    display::print("  System Length(m):              ", xL*L);
    display::print("  Grid Points:                ", ni);
    display::print("  Normalized Time Step:       ", DT);
    display::print("  Total Time Steps:           ", NUM_TS);
    display::print("  Plasma Density:             ", density);
    
    display::print("\n Plasma Characteristics");
    display::print("--------------------------------------------");
    display::print("  Electron Debye Length:      ", LDe);
    display::print("  Ion Debye Length:           ", LDi);
    display::print("  Electron Plasma Frequency:  ", wpe);
    display::print("  Ion Plasma Frequency:       ", wpi);
    display::print("  Electron Thermal Velocity:  ", LDe * wpe, " m/s");
    display::print("  Ion Thermal Velocity:       ", LDi * wpi, " m/s");
    display::print("  Ion Acoustic Velocity:      ", IAW_vel, " m/s");
    
    display::print("\n Normalization");
    display::print("--------------------------------------------");
    display::print("  Normalization Scheme:       ", 
                        normscheme == 1 ? "Electron-scale" :
                        normscheme == 2 ? "Ion-scale" :
                        normscheme == 3 ? "Sub-cycling" :
                        normscheme == 4 ? "Time in electron, space in ion scale" :
                        "Custom");
    
    display::print("  Velocity Normalization:     ",
                        vel_normscheme == 1 ? "Electron thermal velocity" :
                        vel_normscheme == 2 ? "Ion thermal velocity" :
                        "Ion acoustic velocity");
    
    display::print("  Velocity Normalization Factor:  ", vel_norm);
    display::print("  Simulation Time Period:         ", (2 * Const::PI) / W);
    display::print("  Actual Time Step:               ", DT / W);
    display::print("  Total Simulation Time:          ", (DT / W) * NUM_TS, " seconds");
    display::print("  Ion/Electron Thermal Velocity Ratio: ", vel_ratio);
    display::print("  Electron/Ion Thermal Velocity Ratio: ", 1 / vel_ratio);
    
    display::print("\n Solver and Boundary Conditions");
    display::print("--------------------------------------------");
    display::print("  Boundary Condition:  ", bc);
    display::print("  Potential Solver:    ", SolverType);
    display::print("  Solver Tolerance:    ", tolerance);
    display::print("  Max Iterations:      ", max_iteration);
    display::print("  Shape Function:      ", shape);
    
    display::print("\n Diagnostics");
    display::print("--------------------------------------------");
    display::print("  Diagnostics Type:         ", diagtype);
    display::print("  Write Interval:           ", write_interval);
    //display::print("  Write Interval Phase:     ", write_interval_phase);
    display::print("  Sub-cycle Interval:       ", sub_cycle_interval);
   
    
    display::print("\n Collision Properties");
    display::print("--------------------------------------------");
    display::print("  Neutral Gas Density: ", GAS_DENSITY);
    display::print("  Elastic Collision:   ", enable_elastic_collision ? "Enabled" : " Disabled");
    display::print("  Excitation Collision:", enable_excitation_collision ? "Enabled" : " Disabled");
    display::print("  Ionization Collision:", enable_ionization_collision ? "Enabled" : " Disabled");
    display::print(" pion_elastic_collision:", enable_pion_elastic ? "Enabled" : " Disabled");
    display::print(" nion_elastic_collision:", enable_nion_elastic ? "Enabled" : " Disabled");
    display::print("  Detachment Collision:", enable_e_detach_collision ? "Enabled" : " Disabled");
    display::print("  maximum electron collision frequency:", max_electron_coll_freq);
    display::print("  \u03BD * DT :",max_electron_coll_freq*(DT/wpe));


    display::print("\n Magnetic Field");
    display::print("--------------------------------------------");
    display::print("  Bx  : ",  B * sin(theta) * cos(azimuth));
    display::print("  By  : ", B * sin(theta) * sin(azimuth));
    display::print("  Bz  : ", B * cos(theta));

    display::print("\n RF (testing )");
    display::print("--------------------------------------------");
    display::print("  RF cycle  : ", NUM_TS * (DT/W) * fL);
   
    display::print("\n Execution Mode");
    display::print("--------------------------------------------");
    display::print("  Mode: ", num_threads == 1 ? " Serial" : "Parallel (" + std::to_string(num_threads) + " threads)");
    
    int index = 1;
    for(Species &p : species)
    {
        cout << "\n Species (" << index << ") Information ";
        cout << "\n--------------------------------------------";
        cout << "\n  Name:                   " << p.name;
        cout << "\n  Mass:                   " << p.defaultmass;
        cout << "\n  Charge:                 " << p.charge;
        cout << "\n  Temperature:            " << p.temp;
        cout << "\n  Superparticle Weight:   " << p.defaultspwt;
        cout << "\n  Particle Count:         " << p.numparticle;
        cout << "\n  Streaming Velocity:     " << p.vs;
        cout << "\n  Normalized Density:     " << p.fract_den;
        cout << "\n  Posistion Initialization Type:    " << p.initialization_pos;
        cout << "\n  Velocity Initialization Type:    " << p.initialization_vel;
        cout << "\n--------------------------------------------\n";
        index++;
    }

    display::print("test hot cold beam",tau);
    display::print("\n Simulation Ready to Start!\n");
} 

//set normalized parameter.
void Domain:: set_normparam(double LDe, double wpe, double LDi, double wpi)
{
    this->LDe = LDe;
    this->wpe = wpe;
    this->LDi = LDi;
    this->wpi = wpi;  
}

void Domain::set_userdefined_normscheme(double time_scale, double lenght_scale, double energy_scale)
{
    this->time_scale = time_scale;
    this->lenght_scale = lenght_scale;
    this->energy_scale = energy_scale;
}

//new code delete old function if this work
void Domain::set_normscheme()
{
    if (normscheme == 1 || normscheme == 3)
    {
        L = LDe;
        W = wpe;
        if(vel_normscheme == 1)
        {
            vel_norm = LDe*wpe;
        }
        else if(vel_normscheme == 2)
        {
            vel_norm = LDi*wpi;
        }
        else if(vel_normscheme == 3)
        {
            vel_norm = IAW_vel;
        }
        else
        {
            vel_norm = L*W;
        }
        
    }
    else if (normscheme == 2 || normscheme == 4)
    {
        L = LDi;
        W = (normscheme == 2) ? wpi : wpe;  // W remains unchanged for normscheme == 2
        if(vel_normscheme == 1)
        {
            vel_norm = LDe*wpe;
        }
        else if(vel_normscheme == 2)
        {
            vel_norm = LDi*wpi;
        }
        else if(vel_normscheme == 3)
        {
            vel_norm = IAW_vel;
        }
        else
        {
            vel_norm = L*W;
        }

    }
    else if (normscheme == 5)
    {
        LDe = lenght_scale;
        wpe = time_scale;
        L = LDe;
        W = wpe;
        
        vel_norm = L*W;
    
    }
    else
    {
        throw std::invalid_argument("Invalid normscheme flag");
    }
}


//set time.
void Domain::set_time(double DT, int NUM_TS, int write_interval)
{
    this->DT = DT;
    this->NUM_TS = NUM_TS;
    this->write_interval = write_interval;
}

//compute normalized charge density (n_e/N) at grid points.
void Domain::ComputeRho(vector<Species> &species)// "species" is conatiner(vector here) holding instances(many instances) of class "Species"
{
    //":" means range based for loop. here we are iterating over the elements contained in
    // "species" container using a variable sp(adress) and Species is kind of data type(we can also use auto if we we donot know data type )
    // here data type is complex not like int,float char. it is "Species".
    //we use address because we want to amke changes in the original list and donot want to  make a copy.
    //for(Species &sp : species)
    
    // Iterate over each grid point
    for (int i = 0; i < ni; i++) 
    {
        // Initialize density at current grid point to zero
        rho(i) = 0.0;
        // Iterate over each species in the vector
        for (Species &sp : species) 
        {
            // Add contribution of current species to rho at current grid point
            rho(i) += (1 / Const::QE) * (sp.charge * sp.den(i));
        }
    }    
}


//test code

//scatter densities at grid points
void Domain::Scatter(double lc, double value, vec<double> &field)
{
    int i = (int)lc;
    //int i = floor(lc);
    //cout<<ni<<endl;
    double di = lc-i;

    int ngp = round(lc);
    
    if(shape == "CIC")
    {
        field(i) += value*(1-di);
        field(i+1) += value*(di);
    }
    else if(shape == "NGP")
    {
        field(ngp) +=value;
    }
    else if(shape == "TSC")
    {
        if(ngp >0 && ngp < ni-1)
        {
            field(ngp) += value*(0.75 - (lc-ngp)*(lc-ngp));
            field(ngp + 1) += value*( 0.5 * (0.5 + fabs(lc-ngp)) * (0.5 + fabs(lc-ngp)));
            field(ngp - 1) += value*( 0.5 * (0.5 - fabs(lc-ngp)) * (0.5 - fabs(lc-ngp)));
        }
    }
}


//Interpolate density and eletric field value at grid points
// CIC = Cloud in cell
// NGP = Nearest grid point
// TSC = Triangular shape cloud(Quadratic spline)
double Domain::Interpolate(double lc, vec<double> &field)
{
    int i=(int)lc;
    double di = lc-i;
    int ngp = round(lc);
    double val = 0;
    double weight;
    double dis;

    if(shape == "CIC")
    {
        val = field(i)*(1-di) + field(i+1)*(di);
    }
    else if(shape == "NGP")
    {
       val = field(ngp);
    }
    else if(shape == "TSC")
    {
        if(ngp >0 && ngp < ni-1)
        {
            val = field(ngp)*(0.75 - (lc-ngp)*(lc-ngp)) + field(ngp + 1)*( 0.5 * (0.5 + fabs(lc-ngp)) * (0.5 + fabs(lc-ngp))) + field(ngp - 1)*( 0.5 * (0.5 - fabs(lc-ngp)) * (0.5 - fabs(lc-ngp)));
        }
    }
    
    return val;
}

//test code


//this function doesmot needed when the code is normalized
double Domain::XtoL(double pos)
{
    double li = (pos-x0)/dx;
    //cout<<li<<endl;
    return li;
}

double Domain::ComputePE(Species &species)
{
    double pe = 0.0;
    double unnorm = (density * Const::QE * L) / Const::EPS_0;
    vec<double> E_square(ni);
    
    // Calculate E_square = (ef[i] * unnorm)^2 for each i
    for (int i = 0; i < ni; i++)
    {
        E_square(i) = pow(ef(i) * unnorm, 2.0);
    }

    // Perform trapezoidal integration of E_square over L
    double E_int = simpson_38(E_square, L, ni);

    // Calculate potential energy (pe) using the integrated electric field squared
    pe = 0.5 * Const::EPS_0 * E_int;
   
    double Th;
    if(normscheme == 5)
    {
        Th = energy_scale;
    }
    else
    {
        Th = (species.temp*Const::eV)*(species.defaultspwt)*species.numparticle;
    }

    // Calculate total thermal energy (Th) of electrons in the species
    //double Th = (tempE * Const::eV) * species.spwt * species.numparticle;

    pe /= Th;

    return pe;
}


void Domain::filter(vec<double> &field)
{
    vec<double> *filtered_field = &field;
    (*filtered_field)(0) = 0.5*(field(0) + field(1));
    (*filtered_field)(ni-1) = 0.5*(field(ni-2) + field(ni-1));
    for(int i = 1 ; i < ni -1; i++)
    {
        (*filtered_field)(i) = 0.25*(field(i+1) + 2*field(i)+ field(i-1));
    }
}



double Domain::Calculate_alpha(Species &species1, Species &species2)
{
    double n1 = 0.0;
    double n2 = 0.0;
    for (auto &p : species1.part_list)
    {
        n1 += p.spwt;
    }

    for (auto &p : species2.part_list)
    {
        n2 += p.spwt;
    }

    return (n2/n1);

}

double Domain::unirand(double lower_bound, double upper_bound) 
{
    // Create a random number generator
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator

    // Define the distribution range
    std::uniform_real_distribution<> distr(lower_bound, upper_bound);

    // Generate and return a random number in the defined range
    return distr(gen);
}
 
