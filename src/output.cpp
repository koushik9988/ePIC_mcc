#include "output.h"

Output::Output(const std::filesystem::path& outputfolder, Domain& domain) : outputfolder(outputfolder), domain(domain) 
{

    std::filesystem::remove_all(outputfolder);  // Clear previous output
    std::filesystem::create_directories(outputfolder);


    // Convert std::filesystem::path to const char*
    std::string filename = (outputfolder / "result.h5").string();
    const char* c_filename = filename.c_str();
    file = H5File(c_filename, H5F_ACC_TRUNC);


    //file = H5File(outputfolder / "result.h5", H5F_ACC_TRUNC);
    //file = H5File(outputfolder / "result.h5",H5F_ACC_RDWR);
    
    
    if (file.getId() < 0) 
    {
        throw std::runtime_error("Error opening HDF5 file");
    }
    //create groups
    //particle_group1 = file.createGroup("/electron");
    //particle_group2 = file.createGroup("/ion");
    field_data_group = file.createGroup("/fielddata");
    time_group = file.createGroup("/time_var");
    metadata_group = file.createGroup("/metadata");
    metadata_species = file.createGroup("/metadata_species");



    int sp_no = domain.species_no;
    int t_step = int(domain.NUM_TS/domain.write_interval) + 1;
    
    store_ke = Matrix<double>(t_step,sp_no + 2);
    //store_m = Matrix<double>(t_step,sp_no + 1);
    store_m = Matrix<double>(t_step, sp_no * 3 + 1); //new structure to allocate extra memeory for 3-momentum component


}

//-------
void Output::write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT, double density, int save_fig, int normscheme, 
    int subcycleint, double LDe, double LDi, double wpe, double wpi,int spno, double GAS_DENSITY, double max_electron_collision_freq)
{
    // Write metadata attributes within the metadata group
    //Group metadata_group = file.openGroup("/metadata");

    metadata_group.createAttribute("NC", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &NC);
    metadata_group.createAttribute("NUM_TS", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &NUM_TS);
    metadata_group.createAttribute("write_int", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int);
    metadata_group.createAttribute("write_int_phase", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int_phase);
    metadata_group.createAttribute("DT_coeff", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &DT);
    metadata_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &density);
    metadata_group.createAttribute("save_fig", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &save_fig);
    metadata_group.createAttribute("norm_scheme", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &normscheme);
    metadata_group.createAttribute("sub_cycle_interval", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &subcycleint);
    metadata_group.createAttribute("LDe", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &LDe);
    metadata_group.createAttribute("LDi", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &LDi);
    metadata_group.createAttribute("wpe", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &wpe);
    metadata_group.createAttribute("wpi", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &wpi);
    metadata_group.createAttribute("spno", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &spno);
    metadata_group.createAttribute("GAS_DENSITY", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &GAS_DENSITY);
    metadata_group.createAttribute("max_ele_coll_freq", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &max_electron_collision_freq);
    metadata_group.close();
}


void Output::write_species_metadata(std::vector<Species> &species_list)
{
    for (Species& sp : species_list)
    {
        std::string species_group_name = sp.name;
        Group species_group;
   
        species_group = metadata_species.createGroup(species_group_name);
    
        // Write attributes specific to the species
        species_group.createAttribute("name", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.name.c_str());
        species_group.createAttribute("mass", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.mass);
        species_group.createAttribute("charge", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.charge);
        species_group.createAttribute("spwt", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.spwt);
        species_group.createAttribute("temperature", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.temp);
        species_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.fract_den);
        species_group.createAttribute("num_particles", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &sp.numparticle);
        species_group.createAttribute("streaming_velocity", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.vs);

        // Close the species group after writing the metadata
        species_group.close();
    }
}

void Output::write_field_data(int ts)
{
    
    std::string pot_group_name  = "pot";
    std::string efield_group_name  = "efield";

    Group pot_subgroup;
    Group efield_subgroup;

    if (!field_data_group.exists(pot_group_name))
    {
        // Subgroup doesn't exist, create it
        pot_subgroup = field_data_group.createGroup(pot_group_name);
        efield_subgroup = field_data_group.createGroup(efield_group_name);
    }
    else
    {
        // Subgroup already exists, retrieve it
        pot_subgroup = field_data_group.openGroup(pot_group_name);
        efield_subgroup = field_data_group.openGroup(efield_group_name);
    }
    
    //pot_subgroup = field_data_group.createGroup(subgroup_name);
    

    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_pot(Rank, dims_den);
    DataSpace dataspace_efield(Rank, dims_den);

    H5::DataSet dataset_pot = pot_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_pot);
    H5::DataSet dataset_efield = efield_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_efield);

    // Prepare data buffer
    std::vector<double> pot(ni);
    std::vector<double> efield(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        pot[i] = domain.phi(i);
        efield[i] = domain.ef(i);  
    }

    // Write the density data to the dataset
    dataset_pot.write(pot.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_efield.write(efield.data(), H5::PredType::NATIVE_DOUBLE);
}

//-----------------den data----------------------------------------
void Output::write_den_data(int ts,  Species& species)
{
    
    std::string subgroup_name = "den_" + species.name;

    Group den_subgroup;

    // Check if the group already exists in the map
    auto it = den_subgroups.find(species.name);
    if (it != den_subgroups.end()) 
    {
        // Group already exists, retrieve it from the map
        den_subgroup = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        den_subgroup = field_data_group.createGroup(subgroup_name);
        den_subgroups[species.name] = den_subgroup;
    }
    
    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_den(Rank, dims_den);

    H5::DataSet dataset_den = den_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_den);

    // Prepare data buffer
    std::vector<double> density(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        density[i] = species.den(i);
    }

    // Write the density data to the dataset
    dataset_den.write(density.data(), H5::PredType::NATIVE_DOUBLE);
}


// vel data 

void Output::write_vel_data(int ts,  Species& species)
{
    
    std::string subgroup_name = "vel_" + species.name;

    Group vel_subgroup;

    // Check if the group already exists in the map
    auto it = vel_subgroups.find(species.name);
    if (it != vel_subgroups.end()) 
    {
        // Group already exists, retrieve it from the map
        vel_subgroup = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        vel_subgroup = field_data_group.createGroup(subgroup_name);
        vel_subgroups[species.name] = vel_subgroup;
    }
    
    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_vel(Rank, dims_den);

    H5::DataSet dataset_vel = vel_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_vel);

    // Prepare data buffer
    std::vector<double> density(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        density[i] = species.velmesh(i);
    }

    // Write the density data to the dataset
    dataset_vel.write(density.data(), H5::PredType::NATIVE_DOUBLE);
}

//
void Output::write_particle_data(int ts, Species& species)
{
    std::string group_name = "particle_" + species.name;
    Group particle_group;

    // Check if the group already exists in the map
    auto it = particle_groups.find(species.name);
    if (it != particle_groups.end()) 
    {
        // Group already exists, retrieve it from the map
        particle_group = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        particle_group = file.createGroup(group_name);
        particle_groups[species.name] = particle_group;

    }

    // Create datasets for particle positions and velocities
    hsize_t dims_phase[2] = {species.part_list.size(), 4}; //here 2 => 1 for pos and 1 for vel
    //hsize_t dims_vel[2] = {species.part_list.size(), 1};
    hsize_t Rank = 2;

    DataSpace dataspace_phase(Rank, dims_phase);
    //DataSpace dataspace_vel(Rank, dims_vel);

    H5::DataSet dataset_phase = particle_group.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_phase);
    //H5::DataSet dataset_vel = particle_group.createDataSet("vel" + std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_vel);

    // Write particle positions and velocities to the datasets
    std::vector<double> phase_data;
    //std::vector<double> velocities;

    for (Particle& p : species.part_list) 
    {
        phase_data.push_back(p.x);
        //velocity at time t is average of v(t-0.5*dt) and v(t+0.5*dt)
        phase_data.push_back(p.vx);
        phase_data.push_back(p.vy);
        phase_data.push_back(p.vz);
    }

    dataset_phase.write(phase_data.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vel.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
}

void Output::write_ke()
{
    // Define the name for the dataset
    std::string datasetName = "kinetic_energy";

    // Define the dimensions of the dataset
    hsize_t ny = int(domain.NUM_TS/domain.write_interval) + 1; //rows
    hsize_t nx = hsize_t(domain.species_no + 2); //column

    hsize_t dims_energy[2] = {ny, nx}; // Assuming all rows have the same length

    // Create dataspace for the dataset
    DataSpace dataspace_energy(2, dims_energy);

    // Create the dataset within the time_group
    H5::DataSet dataset_energy = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace_energy);

    // Allocate memory
    double* ke_data = new double[nx * ny];

    // Fill the data array
    for (hsize_t i = 0; i < ny; ++i)
    {
        for (hsize_t j = 0; j < nx; ++j)
        {
            ke_data[i * nx + j] = store_ke(i,j);
        }
    } 

    // Write the data to the dataset
    dataset_energy.write(ke_data, H5::PredType::NATIVE_DOUBLE);

    // Deallocate memory after writing
    //delete[] ke_data;
}
//*/
void Output::printmatrix(int row, int col, double **matrix)
{
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            std::cout << matrix[i][j] <<"\t\t";
        }
        std::cout << std::endl;
    }
}

void Output::storeKE_to_matrix(int ts, std::vector<Species> &species_list)
{
    
    auto norm_species = (domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0];

    int k = ts/domain.write_interval;

    store_ke(k,0) = ts * domain.DT;
    //store_ke[][] = domain.ComputePE(norm_species);
    //output.store_ke[int(index/write_interval)][0] = ts * domain.DT;
    //display::print(k);
    int j = 1 ;
    for (Species &sp : species_list)
    {
        store_ke(k,j) = sp.Compute_KE(norm_species);
        j++;
        //cout<<output.store_ke[k][j]<<",";
    }
    store_ke(k,j) = domain.ComputePE(norm_species);
}

//------------
void Output::storem_to_matrix(int ts, std::vector<Species> &species_list)
{
    auto norm_species = (domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0];
    
    int k = ts / domain.write_interval;

    store_m(k, 0) = ts * domain.DT;  // Store time step

    int j = 1;
    for (Species &sp : species_list)
    {
        auto [px, py, pz] = sp.Compute_Momentum(norm_species);

        // Store each component in consecutive columns
        store_m(k, j) = px;
        store_m(k, j + 1) = py;
        store_m(k, j + 2) = pz;
        
        j += 3; // Move to the next set of three columns for the next species
    }
}

void Output::write_m()
{
    // Define the name for the dataset
    std::string datasetName = "momentum";

    // Define the dimensions of the dataset
    hsize_t ny = int(domain.NUM_TS / domain.write_interval) + 1; // Number of rows (time steps)
    hsize_t nx = hsize_t(domain.species_no * 3 + 1); // Number of columns (time + 3 components per species)

    hsize_t dims_momentum[2] = {ny, nx};

    // Create dataspace for the dataset
    DataSpace dataspace_momentum(2, dims_momentum);

    // Create the dataset within the time_group
    H5::DataSet dataset_momentum = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace_momentum);

    // Allocate memory for the data
    double *momentum_data = new double[nx * ny];

    // Fill the data array
    for (hsize_t i = 0; i < ny; ++i)
    {
        for (hsize_t j = 0; j < nx; ++j)
        {
            momentum_data[i * nx + j] = store_m(i, j);
        }
    } 

    // Write the data to the dataset
    dataset_momentum.write(momentum_data, H5::PredType::NATIVE_DOUBLE);

    // Deallocate memory after writing
    //delete[] momentum_data;
}



// Other includes and code

void Output::diagnostics(int ts, std::vector<Species> &species_list)
{
    
    if (domain.diagtype == "off")
    {
        return; // Exit the function if diagnostics are turned off
    }
    
    double max_phi = domain.phi(0);
    double min_phi = max_phi;
    double total_kinetic_energy = 0.0;
    double total_px = 0.0;
    double total_py = 0.0;
    double total_pz = 0.0;
    //double unnorm = Const::K_b*domain.tempE/Const::eV;

    auto norm_species = (domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0];
    double potential_energy;// = domain.ComputePE(norm_species);

    for(int i = 0; i < domain.ni; i++)
    {
        if (domain.phi(i) > max_phi)
        {
            max_phi = domain.phi(i);
        }

        if (domain.phi(i) < min_phi)
        {
            min_phi = domain.phi(i);
        }
    }

    if(domain.SolverType == "direct")
    {
        std::cout << "TS: " << ts << " \t max_phi: " << std::fixed << std::setprecision(2) << (max_phi - domain.phi(0));
    }
    else
    {
        std::cout << "TS: " << ts << "\t" << "norm:" << domain.norm << " delta_phi: " << std::fixed << std::setprecision(2) << (max_phi - domain.phi(0));
    }

    if(domain.bc =="open")
    {
        std::cout<<"\t"<<std::fixed << std::setprecision(precision)<<"rhoL|currentL:"<<domain.vL<<"|"<<domain.I_leftwall<<"\t"<<"rhoR|currentR:"<<domain.vR<<"|"<<domain.I_rightwall;
    }

    if(domain.bc == "open" || domain.enable_ionization_collision  == true)
    {
        for (Species& sp : species_list)
        {
            std::cout << " n_" << std::setw(4) << sp.name << ":" << sp.part_list.size();
        }
    }

    if(domain.normscheme == 5)
    {
        std::cout << std::scientific << std::setprecision(precision);
    }
    else
    {
        std::cout << std::fixed << std::setprecision(precision);
    }
    
    //
    if ((domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full" && Energy_plot == 1)
    {
        potential_energy = domain.ComputePE(norm_species);
        for (Species& sp : species_list)
        {
            double ke = sp.Compute_KE(norm_species);
            total_kinetic_energy += ke;
            std::cout << " KE_" << std::setw(4) << sp.name << ": " << ke;
        }

        std::cout << " Potential energy: " << potential_energy;
        std::cout << " Total_energy: " << total_kinetic_energy + potential_energy;

        for (Species& sp : species_list)
        {
            //total_momentum += sp.Compute_Momentum(norm_species);
            auto [px, py, pz] = sp.Compute_Momentum(norm_species);
            total_px += px;
            total_py += py;
            total_pz += pz;
        }
        std::cout << " px: " << total_px << " py: " << total_py << " pz: " << total_pz<<" P:"<<sqrt(total_px*total_px + total_py*total_py + total_pz*total_pz);
        std::cout<< " delta_g:"<<domain.delta_g;
        
    }
    std::cout << std::endl;

    // Plotting with matplotlibcpp
    static std::vector<double> time;
    static std::vector<double> lenght;
    static std::vector<double> ke;
    static std::vector<double> pe;
    static std::vector<double> total_energie;
    static std::vector<double> phi;
    static std::vector<double> charge_density;
    static std::vector<double> efield;
    static std::vector<double> pot;
    static std::vector<int> num1;
    static std::vector<int> num2;

    //dft plot
    static std::vector<double> k;
    static std::vector<double> rho_k;

    std::vector<double> x;
    std::vector<double> vx;

    std::vector<double> x1;
    std::vector<double> vx1;

    time.push_back(ts * domain.DT);

    if (Energy_plot == 1 && domain.diagtype == "full")
    {
        ke.push_back(total_kinetic_energy);
        pe.push_back(potential_energy);
        total_energie.push_back(total_kinetic_energy + potential_energy);
    }

    if (domain.bc == "open" && domain.diagtype == "full")
    {
        num1.push_back(species_list[0].part_list.size());
        num2.push_back(species_list[1].part_list.size());
    }

    if (phase_plot == 1 && domain.diagtype == "full")
    {
        //x.clear();
        //vx.clear();
        for(auto &part : species_list[species_index].part_list)
        {
            x.push_back(part.x);
            vx.push_back(part.vx);
            //x.push_back(part.vel[0]);
            //vx.push_back(part.vel[1]);
        }
///////////////////////////////////////////////
        //for(auto &part : species_list[2].part_list)
        //{
            //x.push_back(part.x);
            //vx.push_back(part.vx);
            //x.push_back(part.vel[0]);
            //vx.push_back(part.vel[1]);
        //}
    }

    if (Chargedensity_plot == 1 && domain.diagtype == "full")
    {
        lenght.clear();
        charge_density.clear(); 
        for(int i = 0; i <domain.ni ; i++)
        {
            lenght.push_back(i*domain.dx);
            charge_density.push_back(domain.rho(i));
            //charge_density.push_back(species_list[0].den(i));
        }   
    }

    if((Potentialfield_plot == 1 || domain.bc == "open") && domain.diagtype == "full")
    {
        lenght.clear();
        efield.clear();
        pot.clear();
        for(int i = 0; i <domain.ni ; i++)
        {
            lenght.push_back(i*domain.dx);
            efield.push_back(domain.ef(i));
            pot.push_back(domain.phi(i));
        }   
    }

    if(dft_flag == 1 && domain.diagtype == "full")
    {
        k.clear();
        rho_k.clear();
        for(int i = 0; i <domain.ni ; i++)
        {
            k.push_back(domain.dft_k(i));
            rho_k.push_back(domain.dft_value(i));
        }
    }

    plt::ion(); 

    if (Energy_plot == 1 && (domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full")
    {
        plt::figure(1);
        plt::clf();
        if (keflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full")
        {
            plt::named_plot("kinetic energy", time, ke, "r-");
        }
        if (peflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full")
        {
            plt::named_plot("potential energy", time, pe, "b-");
        }
        if (teflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc" ) && domain.diagtype == "full")
        {
            plt::named_plot("total energy", time, total_energie, "g-");
        }
        plt::xlabel("Time");
        plt::ylabel("Energy");
        plt::legend();
    }


    double marker_size = 1.0;
    std::map<std::string, std::string> style_options1 = {{"color", "blue"}, {"marker", "o"}};   
    std::map<std::string, std::string> style_options2 = {{"color", "red"}, {"marker", "o"}};

    if (phase_plot == 1 && domain.diagtype == "full")
    {   
        std::string label1 = species_list[species_index].name ;
        //std::string label2 = species_list[2].name ;

        std::map<std::string, std::string> scatter_keywords1;
        scatter_keywords1["label"] = label1;
        scatter_keywords1["color"] = "black"; // Explicitly set color

        //std::map<std::string, std::string> scatter_keywords2;
        //scatter_keywords2["label"] = label2;
        //scatter_keywords2["color"] = "red"; // Explicitly set color

        // Plot the scatter plots
        plt::figure(2);
        plt::clf();
        plt::scatter(x, vx, marker_size, scatter_keywords1);
        //plt::scatter(x1, vx1, marker_size, scatter_keywords2);
        plt::xlabel("x");
        plt::ylabel("v");

        // Create legend keywords map
        std::map<std::string, std::string> legend_keywords;
        legend_keywords["loc"] = "upper right"; // Location of the legend

        // Apply legend with location
        plt::legend(legend_keywords); 

        // Show the plot
        plt::show(); 
    }


    if (Chargedensity_plot == 1 && domain.diagtype == "full")
    {   
        plt::figure(3);
        plt::clf();
        plt::named_plot("charge-density", lenght, charge_density, "r-");
        plt::xlabel("x");
        plt::ylabel("rho");
        plt::legend(); 
    }

    if((Potentialfield_plot == 1 || domain.bc == "open") && domain.diagtype == "full")
    {
        plt::figure(4);
        plt::clf();
        plt::named_plot("potential", lenght, pot, "r-");
        //plt::named_plot("Electricfield", lenght, efield, "b-");
        plt::xlabel("x");
        //plt::ylim(min_phi,max_phi);
        plt::ylabel("phi/Efield");
        plt::legend(); 
    }

    if(dft_flag == 1 && domain.diagtype == "full")
    {
        plt::figure(6);
        plt::clf();
        plt::named_plot("fourier transformed charge density", k, rho_k, "r-");
        plt::xlabel("k");
        plt::ylabel("rho_k");
        plt::legend(); 
    }
    
    if(domain.bc == "open" && domain.diagtype == "full")
    {   
        plt::figure(5);
        plt::clf();
        plt::named_plot("electron", time, num1, "r-");
        plt::named_plot("ion", time, num2, "b-");
        //plt::named_plot("Electricfield", lenght, efield, "b-");
        plt::xlabel("time");
        plt::ylabel("simulation particle");
        plt::legend();

    }

    plt::pause(0.1);

    // Show plot
    plt::show();
}
