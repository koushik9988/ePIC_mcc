/* */
#include "iniparser.h"
#include <iostream>
#include "field.h"
#include "init.h"
#include "domain.h"
#include "species.h"
#include <cmath>
#include <fstream>
#include <time.h>
#include <chrono>
#include "output.h"
#include "extrafun.h"
#include <thread>
#include <string>
#include "collision.h"
#include <tuple>

using namespace std; 
using namespace display;

int main( int argc , char *argv[]) 
{     
    if(argc<2)
    {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }

    
    //parsing input.ini file and storing values
    const std::string filename = argv[1];
    
    int num_threads;
    if (argc>2)
    {
        num_threads = atoi(argv[2]);
    }
    else
    {
        num_threads = std::thread::hardware_concurrency();
    }

    if(num_threads > std::thread::hardware_concurrency())
    {
        num_threads = std::thread::hardware_concurrency();
    }
    if(num_threads == 1)
    {
       //display::print("runing serial code");
    }

    //display::print("running with ",num_threads," threads");
    
    auto iniData = INIParser::parse(filename);

    auto species_section = iniData["species"];
    int species_no = species_section.size();

    //output folder
    std::string outputfolder = INIParser::getString(iniData["file"],"output");

    //grid/domain
    int ni = INIParser::getInt(iniData["domain"], "ni");//grid no
    //int ni = NC+1; // no of grid points is one more than cell no
    double x0 = INIParser::getInt(iniData["domain"],"x0");

    //diagnostic
    int save_fig = INIParser::getDouble(iniData["diagnostics"],"save_fig");
    int write_interval = INIParser::getInt(iniData["diagnostics"],"write_interval");
    int write_interval_phase = INIParser::getInt(iniData["diagnostics"],"write_interval_phase");
    int write_diagnostics = INIParser::getInt(iniData["diagnostics"],"write_diagnostics");
    int sub_cycle_interval = INIParser::getInt(iniData["diagnostics"],"sub_cycle_interval");
    int precision = INIParser::getInt(iniData["diagnostics"],"precision");
    int write_flag = INIParser::getInt(iniData["diagnostics"],"write_flag");
    std::string diagtype = INIParser::getString(iniData["diagnostics"],"diagtype");

    //time
    int NUM_TS = INIParser::getInt(iniData["time"], "NUM_TS");
    double DT_coeff = INIParser::getDouble(iniData["time"],"DT_coeff");
    
    //simulation parameter
    int see_rate = INIParser::getInt(iniData["simulation"],"see_rate");
    std:: string bc = INIParser::getString(iniData["simulation"],"bc");
    std:: string shapefunction = INIParser::getString(iniData["simulation"],"shapefunction");
    std::string push_parallal = INIParser::getString(iniData["simulation"],"push_parallal");
    std::string deposit_parallal = INIParser::getString(iniData["simulation"],"deposit_parallal");
    double den = INIParser::getDouble(iniData["simulation"],"density");
    double tempwall = INIParser::getInt(iniData["simulation"],"tempwall");
    int ionfixed = INIParser::getInt(iniData["simulation"],"ionfixed");


    //collision
    std::string elastic_flag = INIParser::getString(iniData["collision"],"elastic");
    std::string excitation_flag = INIParser::getString(iniData["collision"],"excitation");
    std::string ionization_flag = INIParser::getString(iniData["collision"],"ionization");
    double GAS_DENSITY = INIParser::getDouble(iniData["collision"],"GAS_DENSITY");
    std::string collgroup = INIParser::getString(iniData["collision"],"collgroup");
    auto collisionPairs = INIParser::parseCollGroup(collgroup);

    
    //normalization
    int norm_scheme = INIParser::getInt(iniData["normalization"],"norm_scheme");
    int vel_normscheme = INIParser::getInt(iniData["normalization"],"vel_norm_scheme");
    double L_scale = INIParser::getDouble(iniData["normalization"],"lenght_scale");
    std::string T_scale = INIParser::getString(iniData["normalization"],"time_scale");

    //potential solver
    double tolerance = INIParser::getDouble(iniData["solver"],"tolerance");
    double max_iteration = INIParser::getInt(iniData["solver"],"max_iteration");
    std::string SolverType = INIParser::getString(iniData["solver"],"solvertype");
    
    //visual plot flag(using matplotlibcpp.h) 
    int Energyplot_flag = INIParser::getInt(iniData["visualplot"],"Energy_plot");
    int chargeplot_flag = INIParser::getInt(iniData["visualplot"],"Chargedensity_plot");
    int Potentialfieldplot_flag = INIParser::getInt(iniData["visualplot"],"Potentialfield_plot");
    int keflag = INIParser::getInt(iniData["visualplot"],"keflag");
    int peflag = INIParser::getInt(iniData["visualplot"],"peflag");
    int teflag = INIParser::getInt(iniData["visualplot"],"teflag");
    int phaseplot_flag = INIParser::getInt(iniData["visualplot"],"phase_plot");
    int species_index = INIParser::getInt(iniData["visualplot"],"species_index");
    int dft_flag = INIParser::getInt(iniData["visualplot"],"dft_rho");

    //vector to store species data
    std::vector<std::string> names;
    std::vector<double> mass;
    std::vector<int> nParticles;
    std::vector<double> temps;
    std::vector<int> charge_signs;
    std::vector<double> frac_densities;
    std::vector<double> normden;
    std::vector<double> vs;
    std::vector<std::string> pos_init;

    names.reserve(species_no);
    mass.reserve(species_no);
    nParticles.reserve(species_no);
    temps.reserve(species_no);
    charge_signs.reserve(species_no);
    frac_densities.reserve(species_no);
    normden.reserve(species_no);
    vs.reserve(species_no);
    pos_init.reserve(species_no);

    // First pass: Parse the species section and populate vectors
    for (const auto& species_entry : species_section)
    {
        const std::string& line = species_entry.second;
        std::vector<std::string> tokens = INIParser::split(line, ',');

        if (tokens.size() == 8)
        {  
            names.push_back(tokens[0]);                      // Species name
            mass.push_back(std::stod(tokens[1]));
            nParticles.push_back(std::stoi(tokens[2]));      // Number of particles
            temps.push_back(std::stod(tokens[3]));           // Temperature
            charge_signs.push_back(std::stoi(tokens[4]));    // Charge sign (integer -1 or 1)
            frac_densities.push_back(std::stod(tokens[5]));  // Fractional density
            vs.push_back(std::stod(tokens[6]));  // Fractional density
            pos_init.push_back(tokens[7]);
        }
    }

    double k = 0;
    for(int i = 0 ; i < species_no; i++)
    {
        k += (-charge_signs[i])*frac_densities[i];
    }

    //display::print(k);

    normden[1] = den;
    normden[0] = den/k;
    
    for(int i = 2 ;i < species_no; i++)
    {
        normden[i] = frac_densities[i]*normden[0];
    }

    //noramlizing quantity(electron)
    double LDe = sqrt((Const::EPS_0*Const::K_b*temps[0]*Const::EV_to_K)/(normden[0]*Const::QE*Const::QE)); // Electron Debye Length   
    double wpe = sqrt((normden[0]*Const::QE*Const::QE)/(mass[0]*Const::EPS_0)); // Total Electron Plasma Frequency
    double wpi = sqrt((normden[1]*Const::QE*Const::QE)/(mass[1]*Const::EPS_0)); //ion timescale
    double LDi = sqrt((Const::EPS_0*Const::K_b*temps[1]*Const::EV_to_K)/(normden[1]*Const::QE*Const::QE)); // ion Debye Length
    double CS = sqrt(temps[0]*Const::K_b*Const::EV_to_K/mass[1]); // Ion acoustic speed

    double vthe = LDe*wpe;
    double vthi = LDi*wpi;

    //user defined scales
    double energy_scale = 1;
    double time_scale;
    double lenght_scale = L_scale;
    if(T_scale == "omegape")
    {
        time_scale = wpe; 
    }
    else if(T_scale == "omegapi")
    {
        time_scale = wpi; 
    }
    else if(T_scale != "omegape" || T_scale != "omegapi")
    {
        time_scale = std::stod(T_scale);
    }
 
    double dx;
    double DT;
    double stepSize;
    if(norm_scheme == 2 || norm_scheme == 4)
    {
        stepSize = LDi;
        dx = stepSize/LDi;
        DT = DT_coeff*(1.0/wpi);
        DT = wpi*DT;
    }
    if(norm_scheme == 1 || norm_scheme == 3)
    {
        stepSize = LDe;
        dx = stepSize/LDe;
        DT = DT_coeff*(1.0/wpe);
        DT = wpe*DT;
    }
    
    if(norm_scheme == 5)
    {
        stepSize = lenght_scale;
        dx = stepSize/lenght_scale;
        DT = DT_coeff*(1.0/time_scale);
        DT = time_scale*DT;
    }
    

    Domain domain(x0,dx,ni);

    domain.bc  = bc;
    domain.vel_ratio = vthi/vthe;
    domain.see_rate = see_rate;
    domain.tempwall = tempwall;
    domain.shape = shapefunction;
    domain.SolverType = SolverType;
    domain.tolerance = tolerance;
    domain.num_threads = num_threads;
    domain.push_parallal = string_to_bool(push_parallal);
    domain.deposit_parallal = string_to_bool(deposit_parallal);
    domain.species_no = species_no;
    domain.density = den;
    domain.IAW_vel = CS;
    if(domain.num_threads == 1)
    {
        domain.push_parallal = false;
        domain.deposit_parallal = false;
    }
    domain.max_iteration = max_iteration;
    domain.normscheme = norm_scheme;
    domain.vel_normscheme = vel_normscheme;
    domain.diagtype = diagtype;
    domain.sub_cycle_interval = sub_cycle_interval;
    domain.set_normparam(LDe,wpe,LDi,wpi);
    domain.set_userdefined_normscheme(time_scale,lenght_scale,energy_scale);
    domain.set_time(DT,NUM_TS,write_interval);
    domain.set_normscheme();
    domain.vL =  0;
    domain.vR = 0;
    domain.I_leftwall = 0;
    domain.I_rightwall = 0;
    //collision
    domain.GAS_DENSITY = GAS_DENSITY;
    domain.enable_elastic_collision = string_to_bool(elastic_flag);    // Flag for elastic collisions
    domain.enable_excitation_collision = string_to_bool(excitation_flag); // Flag for excitation collisions
    domain.enable_ionization_collision = string_to_bool(ionization_flag); // Flag for ionization collisions
    domain.delta_g = 0;

    domain.ionfixed = ionfixed;
    
    std::vector<Species> species_list;
 
    for (int i = 0 ;i < species_no; i++)
    {
        double computed_spwt = 0;

        computed_spwt = normden[i] * domain.xL * domain.L / nParticles[i];
        //print(domain.L);

        species_list.emplace_back(names[i], mass[i], charge_signs[i]*Const::QE, computed_spwt, temps[i], nParticles[i],vs[i],frac_densities[i], pos_init[i], domain);
    }
 
    
    CollisionHandler ElectronNeutralCollision(domain);
    //collision
    //pre-calculate electron cross-section for energy level(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    ElectronNeutralCollision.set_electron_cross_sections();
    //Calculate total cross-section for energy levels(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    ElectronNeutralCollision.calc_total_cross_sections();
    
    domain.max_electron_coll_freq = ElectronNeutralCollision.max_electron_coll_freq();
    domain.display(species_list);
    print("Press Enter to continue...");
    std::cin.get();

    Output output(outputfolder,domain);

    output.precision = precision;
    output.Energy_plot =  Energyplot_flag ;
    output.Potentialfield_plot = Potentialfieldplot_flag;
    output.Chargedensity_plot = chargeplot_flag ;
    output.keflag = keflag;
    output.peflag = peflag;
    output.teflag = teflag;
    output.phase_plot = phaseplot_flag;
    output.dft_flag = dft_flag;
    output.species_index = species_index;
    output.write_metadata(ni,NUM_TS,write_interval,write_interval_phase,DT_coeff,den,save_fig,domain.normscheme,
        domain.sub_cycle_interval,LDe,LDi,wpe,wpi,species_no,GAS_DENSITY, domain.max_electron_coll_freq);
    output.write_species_metadata(species_list);

    auto start_time = std::chrono::high_resolution_clock::now();

    //initializing the species by creating instances of Init class
    for(Species &sp : species_list)
    {
        Init init(sp,domain);
    } 
    
    //FieldSolve class instance declaration
    FieldSolve fieldsolver(domain);
    
    //scatter species to mesh/grid 
    for (Species &sp:species_list)
	{
		sp.ScatterSpecies();
	}

    //compute charge density on the grid
    domain.ComputeRho(species_list);

    //initial potential and field calculation
    fieldsolver.PotentialSolver(0);
    fieldsolver.CalculateEfield();

    //rewind species velocity by half a time step
    for (Species &sp:species_list)
	{
        sp.Rewind_species();
    }


    //--------------MAIN LOOP----------------------- 
    for(int ts = 0 ; ts < NUM_TS + 1; ts++)
    {
        
        //domain.vL =  5*(Const::eV/(Const::K_b*Const::EV_to_K));
        if(ts%5 == 0)
        {
            //GenerateParticlePairs(species_list[0],species_list[1],1,domain);
        }
        
        for (Species &sp:species_list)
		{
			sp.ScatterSpecies();
            sp.ScatterVel_serial();
		}


        domain.ComputeRho(species_list);
        
        fieldsolver.PotentialSolver(ts);
        fieldsolver.CalculateEfield();

        //---------particle-mover-------------------------------------------- 
        for (Species &sp:species_list)
		{
            
            if (sp.name == "ion" && domain.ionfixed == 1) continue;
            if(domain.normscheme == 1 || domain.normscheme == 2 || domain.normscheme == 4 || domain.normscheme == 5)
            {
                sp.Push_species(species_list,1);
            }
            if(domain.normscheme == 3)
            {
                if(sp.name == "electron")
                {
                    sp.Push_species(species_list,1);
                }
                else
                {
                    if(ts%domain.sub_cycle_interval == 0)
                    {
                        sp.Push_species(species_list,sub_cycle_interval);
                    }
                }
            }
        }

        if(domain.enable_elastic_collision || domain.enable_excitation_collision || domain.enable_ionization_collision)
        {
            for (const auto& [first, second] : collisionPairs)
            {
                ElectronNeutralCollision.handle_collisions(species_list[first], species_list[second]);  
            }
        }
        
        if(ts%write_interval== 0)
        {
            if(write_flag == 1 || write_flag == 2)
            {    
                output.write_field_data(ts);
                output.storeKE_to_matrix(ts,species_list);
                output.storem_to_matrix(ts,species_list);
                for(Species &sp:species_list)
                {
                    output.write_den_data(ts,sp);
                    output.write_vel_data(ts,sp);
                }    
            }
        }
        
        if(ts%write_interval_phase == 0)
        {
            if(write_flag == 1 || write_flag == 3)
            {
                for(Species &sp:species_list)
                {
                    output.write_particle_data(ts,sp);
                }
            }
        }
        
        if(ts%write_diagnostics == 0)
        {
            output.diagnostics(ts,species_list);
        }
    }
    
    output.write_ke();
    output.write_m();

    if(domain.diagtype != "off")
    {
        print("average no of electron crossed one or more than one cell per time step : ",domain.ele_cross/NUM_TS);
        print("average no of ion crossed one or more than  cell per time step : ",domain.ion_cross/NUM_TS,"\n");
        if(domain.ele_cross > 0 )
        {print("average no of cell crossed by electron: ",domain.crossed_cellno_ele/domain.ele_cross,"\n");}
        if(domain.ion_cross > 0 )
        {print("average no of cell crossed by ion: ",domain.crossed_cellno_ion/domain.ion_cross,"\n");}
    
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    if(domain.diagtype != "off")
    {
        std::cout << "Elapsed time: " << elapsed_time.count() << "seconds." <<"or "<< elapsed_time.count()/60<<"minutes"<<std::endl;
    }
    

    print("Press Enter to close the simulation...");
    std::cin.get();

    return 0;
}
