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
#include "emitter.h"

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
    double voltage_left = INIParser::getDouble(iniData["simulation"],"voltage_left");
    double voltage_right = INIParser::getDouble(iniData["simulation"],"voltage_right");



    //collision
    std::string elastic_flag = INIParser::getString(iniData["collision"],"elastic");
    std::string excitation_flag = INIParser::getString(iniData["collision"],"excitation");
    std::string ionization_flag = INIParser::getString(iniData["collision"],"ionization");
    std:: string pion_elastic_flag = INIParser::getString(iniData["collision"],"pion_elastic");
    std:: string e_detach_collision_flag = INIParser::getString(iniData["collision"],"e_detach_collision");
    double GAS_DENSITY = INIParser::getDouble(iniData["collision"],"GAS_DENSITY");
    std::string collgroup = INIParser::getString(iniData["collision"],"collgroup");
    auto collisionPairs = INIParser::parseCollGroup(collgroup);
    std::string GAS_TYPE = INIParser::getString(iniData["collision"],"GAS_TYPE");

    
    //normalization
    int norm_scheme = INIParser::getInt(iniData["normalization"],"norm_scheme");
    int vel_normscheme = INIParser::getInt(iniData["normalization"],"vel_norm_scheme");
    double L_scale = INIParser::getDouble(iniData["normalization"],"lenght_scale");
    std::string T_scale = INIParser::getString(iniData["normalization"],"time_scale");

    //potential solver
    double tolerance = INIParser::getDouble(iniData["solver"],"tolerance");
    double max_iteration = INIParser::getInt(iniData["solver"],"max_iteration");
    std::string SolverType = INIParser::getString(iniData["solver"],"solvertype");


    //external field
    double B = INIParser::getDouble(iniData["ExternalField"],"B");
    double theta = (M_PI / 180.0) * INIParser::getDouble(iniData["ExternalField"],"theta");
    double azimuth = (M_PI / 180.0) * INIParser::getDouble(iniData["ExternalField"],"azimuth");
    
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
    int ke_components = INIParser::getInt(iniData["visualplot"],"ke_components");
    int coll_freq_plot = INIParser::getInt(iniData["visualplot"],"coll_freq_plot");

    
    //parse species data
    const auto& SpeciesSection = iniData.at("Species");

    // Get the number of Species
    int species_no = INIParser::getInt(SpeciesSection, "count");

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
    std::vector<std::string> vel_init;

    names.reserve(species_no);
    mass.reserve(species_no);
    nParticles.reserve(species_no);
    temps.reserve(species_no);
    charge_signs.reserve(species_no);
    frac_densities.reserve(species_no);
    normden.reserve(species_no);
    vs.reserve(species_no);
    pos_init.reserve(species_no);
    vel_init.reserve(species_no);

    std::vector<SpeciesParams> param_list;
    param_list.reserve(species_no);

    for (int i = 0; i < species_no; ++i)
    {
        std::string prefix = "species_" + std::to_string(i) + ".";

        SpeciesParams params;
        params.name = INIParser::getString(SpeciesSection, prefix + "name");
        params.mass = INIParser::getDouble(SpeciesSection, prefix + "mass");
        params.num = INIParser::getInt(SpeciesSection, prefix + "num");
        params.temp = INIParser::getDouble(SpeciesSection, prefix + "temp");
        params.charge_sign = INIParser::getInt(SpeciesSection, prefix + "charge_sign");
        params.normden = INIParser::getDouble(SpeciesSection, prefix + "normden");
        params.vs = INIParser::getDouble(SpeciesSection, prefix + "vs");
        params.loadtype_pos = INIParser::getString(SpeciesSection, prefix + "loadtype_pos");
        params.loadtype_vel = INIParser::getString(SpeciesSection, prefix + "loadtype_vel");
        param_list.push_back(params);
    }

   for (const auto& p : param_list)
   {
        names.push_back(p.name);                      // Species name
        mass.push_back(p.mass);
        nParticles.push_back(p.num);      // Number of particles
        temps.push_back(p.temp);           // Temperature
        charge_signs.push_back(p.charge_sign);    // Charge sign (integer -1 or 1)
        frac_densities.push_back(p.normden);  // Fractional density
        vs.push_back(p.vs);  // Fractional density
        pos_init.push_back(p.loadtype_pos);
        vel_init.push_back(p.loadtype_vel);
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
    domain.vL = voltage_left;
    domain.vR = voltage_right;
    domain.I_leftwall = 0;
    domain.I_rightwall = 0;
    domain.phi_norm = Const::QE/(Const::K_b*Const::EV_to_K*temps[0]); // Normalization factor for potential
    //collision
    domain.GAS_DENSITY = GAS_DENSITY;
    domain.enable_elastic_collision = string_to_bool(elastic_flag);    // Flag for elastic collisions
    domain.enable_excitation_collision = string_to_bool(excitation_flag); // Flag for excitation collisions
    domain.enable_ionization_collision = string_to_bool(ionization_flag); // Flag for ionization collisions
    domain.enable_pion_elastic = string_to_bool(pion_elastic_flag);
    domain.enable_e_detach_collision = string_to_bool(e_detach_collision_flag);
    domain.GAS_TYPE = GAS_TYPE;
    domain.delta_g = 0;
    domain.N_ecoll = 0;
    domain.N_ioncoll = 0;
    domain.N_negioncoll = 0;

    domain.electronegativity = 0;
    
    domain.ionfixed = ionfixed;


    domain.B = B;
    domain.theta = theta;
    domain.azimuth = azimuth;
    
    std::vector<Species> species_list;
 
    for (int i = 0 ;i < species_no; i++)
    {
        double computed_spwt = 0;

        computed_spwt = normden[i] * domain.xL * domain.L / nParticles[i];
        //print(domain.L);

        species_list.emplace_back(names[i], mass[i], charge_signs[i]*Const::QE, computed_spwt, temps[i], nParticles[i],vs[i],frac_densities[i], pos_init[i], vel_init[i], domain);
    }

    ///
    //double beam_temp = 0.5 * species_list[2].mass * pow(species_list[2].vs * domain.vel_norm, 2) / Const::QE;

    //domain.tau = (sqrt(normden[0]*species_list[2].mass))/(sqrt(normden[2]*species_list[0].mass)) * (species_list[2].temp/beam_temp);
    ///
    
    CollisionHandler *ElectronNeutralCollision = nullptr;

    if(domain.GAS_TYPE == "H")
    {
        if(domain.enable_pion_elastic)
        {
            ElectronNeutralCollision = new CollisionHandler(domain, HYDROGEN);
        }
        else
        {
            ElectronNeutralCollision = new CollisionHandler(domain, HYDROGEN, "../cross_section");
        }
        
    }
    else if(domain.GAS_TYPE == "AR")
    {
        ElectronNeutralCollision = new CollisionHandler(domain, ARGON, "../cross_section");
    }

    
    //collision
    //pre-calculate electron cross-section for energy level(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    ElectronNeutralCollision->set_electron_cross_sections();
    ElectronNeutralCollision->set_pion_cross_sections();
    //Calculate total cross-section for energy levels(DE_CS*1,DE_CS*2,DE_CS*3........DE_CS*(CS_RANGES-1))
    ElectronNeutralCollision->calc_total_cross_sections();
    
    domain.max_electron_coll_freq = ElectronNeutralCollision->max_electron_coll_freq();

    

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
    output.ke_components = ke_components;
    output.coll_freq_plot = coll_freq_plot;
    output.write_metadata(ni,NUM_TS,write_interval,write_interval_phase,DT_coeff,den,save_fig,domain.normscheme,
        domain.sub_cycle_interval,LDe,LDi,wpe,wpi,species_no,GAS_DENSITY, domain.max_electron_coll_freq);
    output.write_species_metadata(species_list);


    //-----emitter--
    std::vector<Emitter> emitters;
    const auto& emitterSection = iniData.at("Emitters");

    // Get the number of emitters
    int emitter_count = INIParser::getInt(emitterSection, "count");
    emitters.reserve(emitter_count);

    for (int i = 0; i < emitter_count; i++)
    {
        std::string prefix = "emitter_" + std::to_string(i) + ".";

        EmitterParams params;
        params.emitter_loc = INIParser::getDouble(emitterSection, prefix + "emitter_loc");
        params.temp = INIParser::getDouble(emitterSection, prefix + "temp");
        params.numparticle = INIParser::getInt(emitterSection, prefix + "numparticle");
        params.vdrift = INIParser::getDouble(emitterSection, prefix + "drift_velocity");
        params.species_id1 = INIParser::getInt(emitterSection, prefix + "species_id1");
        params.species_id2 = INIParser::getInt(emitterSection, prefix + "species_id2");
        emitters.emplace_back(params,domain);
    }
    ///---emitter----

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

    //double nu = ElectronNeutralCollision->average_collision_frequency(species_list[0]);
    //display::print("avaerage collision frequecny : ", nu);

    //--------------MAIN LOOP----------------------- 
    for(int ts = 0 ; ts < NUM_TS + 1; ts++)
    {
        
        
        domain.electronegativity = domain.Calculate_alpha(species_list[0], species_list[2]);
        domain.avg_coll_freq = ElectronNeutralCollision->average_collision_frequency(species_list[1]);

        

        //domain.collision_rate.display();
        
        int start = 0;
        int end = NUM_TS;
        if (ts >= start && ts < end)
        {
            for (auto &emit : emitters)
            {
                emit.inject(species_list);
            }
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

        if(domain.enable_elastic_collision || domain.enable_excitation_collision || domain.enable_ionization_collision || domain.enable_pion_elastic || domain.enable_e_detach_collision)
        {
            for (const auto& [first, second] : collisionPairs)
            {
                ElectronNeutralCollision->handle_collisions(species_list[first], species_list[second],species_list[0]);
                ElectronNeutralCollision->coll_rate(species_list[first], species_list[second]);  
            }
            
        }

        //------------------------------------------------------------------------
        
        
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
                    output.write_collrate_data(ts,sp);
                    output.write_vel_data(ts,sp);

                }    
            }
            output.write_avg_collision_freq(ts);
            output.write_alpha_vs_time(ts);
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
    
    //print(domain.N_pcoll);
    //print("total sim time:",NUM_TS*(DT_coeff/wpe));
    double z1 = (domain.N_ecoll/(nParticles[0]*NUM_TS*(DT_coeff/domain.W)))/domain.W;
    double z2 = (domain.N_ioncoll/(nParticles[1]*NUM_TS*(DT_coeff/domain.W)))/domain.W;
    double z3 = (domain.N_negioncoll/((nParticles[2] + nParticles[3])*NUM_TS*(DT_coeff/domain.W)))/domain.W;

    display::print("Average electron collision frequency: ", z1);
    display::print("Average ion collision frequency: ", z2);
    display::print("Average negative ion collision frequency: ", z3);

    output.write_extra(z1,z2,z3);

    print("Press Enter to close the simulation...");
    std::cin.get();

    delete ElectronNeutralCollision;

    return 0;
}
