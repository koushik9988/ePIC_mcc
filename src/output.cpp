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
    
    store_ke = Matrix<double>(t_step,3*sp_no + 2);
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
    std::vector<std::string> species_order;
    for (Species& sp : species_list)
    {
        std::string species_group_name = sp.name;
        Group species_group;
   
        species_group = metadata_species.createGroup(species_group_name);
    
        // Write attributes specific to the species
        //species_group.createAttribute("name", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.name.c_str());
        species_group.createAttribute("mass", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.defaultmass);
        species_group.createAttribute("charge", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.charge);
        species_group.createAttribute("spwt", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.defaultspwt);
        species_group.createAttribute("temperature", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.temp);
        species_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.fract_den);
        species_group.createAttribute("num_particles", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &sp.numparticle);
        species_group.createAttribute("streaming_velocity", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &sp.vs);
        //species_group.createAttribute("pos_init", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.initialization_pos.c_str());
        //species_group.createAttribute("vel_init", PredType::C_S1, DataSpace(H5S_SCALAR)).write(PredType::C_S1, sp.initialization_vel.c_str());
        StrType strdatatype(PredType::C_S1, H5T_VARIABLE); // Define variable-length string type

        species_group.createAttribute("name", strdatatype, DataSpace(H5S_SCALAR)).write(strdatatype, sp.name);
        species_group.createAttribute("pos_init", strdatatype, DataSpace(H5S_SCALAR)).write(strdatatype, sp.initialization_pos);
        species_group.createAttribute("vel_init", strdatatype, DataSpace(H5S_SCALAR)).write(strdatatype, sp.initialization_vel);


        // Add species name to species_order vector
        species_order.push_back(sp.name);
        // Close the species group after writing the metadata
        species_group.close();
    }
    // Create an attribute to store the species order as in species_list, use hsize_t array for the dimension
    hsize_t dims[1] = {species_order.size()};
    DataSpace order_space(H5S_SIMPLE, dims);

    // Define a variable-length string type for the species names as each species name have diffrent lenght strings
    StrType str_type(PredType::C_S1, H5T_VARIABLE);

    Attribute order_attr = metadata_species.createAttribute("species_order", str_type, order_space);

    std::vector<const char*> species_name_pointers;
    for (const auto& name : species_order)
    {
        species_name_pointers.push_back(name.c_str());
    }

    // Write the array of species names as a string array (variable-length strings)
    order_attr.write(str_type, species_name_pointers.data());
}

void Output::write_field_data(int ts)
{
    
    std::string pot_group_name  = "pot";
    std::string efield_group_name  = "efield";
    //std::string coll_rate_group_name = "coll_rate";

    Group pot_subgroup;
    Group efield_subgroup;
    //Group collrate_subgroup;

    if (!field_data_group.exists(pot_group_name))
    {
        // Subgroup doesn't exist, create it
        pot_subgroup = field_data_group.createGroup(pot_group_name);
        efield_subgroup = field_data_group.createGroup(efield_group_name);
        //collrate_subgroup = field_data_group.createGroup(coll_rate_group_name);
    }
    else
    {
        // Subgroup already exists, retrieve it
        pot_subgroup = field_data_group.openGroup(pot_group_name);
        efield_subgroup = field_data_group.openGroup(efield_group_name);
        //collrate_subgroup = field_data_group.openGroup(coll_rate_group_name);
    }
    
    //pot_subgroup = field_data_group.createGroup(subgroup_name);
    

    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_pot(Rank, dims_den);
    DataSpace dataspace_efield(Rank, dims_den);
    //DataSpace dataspace_collrate(Rank, dims_den);

    H5::DataSet dataset_pot = pot_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_pot);
    H5::DataSet dataset_efield = efield_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_efield);
    //H5::DataSet dataset_collrate = collrate_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_collrate);

    // Prepare data buffer
    std::vector<double> pot(ni);
    std::vector<double> efield(ni);
    //std::vector<double> collrate(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        pot[i] = domain.phi(i);
        efield[i] = domain.ef(i);
        //collrate[i] = domain.collision_rate(i);  
    }

    // Write the density data to the dataset
    dataset_pot.write(pot.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_efield.write(efield.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_collrate.write(collrate.data(), H5::PredType::NATIVE_DOUBLE);
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


//##############################
void Output::write_collrate_data(int ts, Species& species)
{
    // Subgroup name: coll_rate_<species>
    std::string subgroup_name = "coll_rate_" + species.name;

    Group collrate_subgroup;

    // Check if the group already exists in the map
    auto it = collrate_subgroups.find(species.name);
    if (it != collrate_subgroups.end()) 
    {
        // Group already exists, retrieve it from the map
        collrate_subgroup = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        collrate_subgroup = field_data_group.createGroup(subgroup_name);
        collrate_subgroups[species.name] = collrate_subgroup;
    }

    hsize_t ni = domain.ni;
    hsize_t dims[1] = {ni};
    hsize_t Rank = 1;

    DataSpace dataspace(Rank, dims);

    // Create dataset for this timestep
    H5::DataSet dataset_collrate = collrate_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace);

    // Prepare buffer
    std::vector<double> collrate(ni);

    // Fill with values from domain
    for (int i = 0; i < ni; ++i) 
    {
        collrate[i] = species.coll_rate(i);
    }

    // Write to HDF5
    dataset_collrate.write(collrate.data(), H5::PredType::NATIVE_DOUBLE);
}

//##############################


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
    hsize_t dims_phase[2] = {species.part_list.size(), 5}; //here 4 => 1 for pos and 3 for vel + 1 extra for particle id
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
        phase_data.push_back(p.id); // Store particle id
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
    hsize_t nx = hsize_t(3*domain.species_no + 2); //column

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
    // Determine normalization species
    auto norm_species = (domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0];

    // Timestep index for storage
    int k = ts / domain.write_interval;

    // Store timestep value
    store_ke(k, 0) = ts * domain.DT;

    //Store kinetic energy for each species (3 components each)
    int j = 1;  // Start after timestep column
    for (Species &sp : species_list)
    {
        vec<double> ke = sp.Compute_KE(norm_species);  // Returns KE vector [x, y, z]
        store_ke(k, j)     = ke(0);  // KE in x-direction
        store_ke(k, j + 1) = ke(1);  // KE in y-direction
        store_ke(k, j + 2) = ke(2);  // KE in z-direction
        j += 3;  // Increment by 3 for next species
    }

    // Store potential energy in the last two column
    store_ke(k, j) = domain.ComputePE(norm_species);
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

void Output::write_avg_collision_freq(int ts)
{
    std::string datasetName = "avg_collision_freq";
    
    hsize_t ny = int(domain.NUM_TS / domain.write_interval) + 1; // rows (time steps)
    hsize_t dims[2] = {ny, 2}; // 2 columns: [time, frequency]

    // Create dataspace for full dataset
    DataSpace dataspace(2, dims);

    // Create or open the dataset
    H5::DataSet dataset;
    if (!time_group.exists(datasetName))
    {
        dataset = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    }
    else
    {
        dataset = time_group.openDataSet(datasetName);
    }

    // Prepare the data for this timestep
    int k = ts / domain.write_interval;
    double data[1][2];
    data[0][0] = ts * domain.DT;                    // Time
    data[0][1] = domain.avg_coll_freq / domain.W;   // Average collision frequency

    // Define hyperslab in the dataset where this row will be written
    hsize_t offset[2] = {static_cast<hsize_t>(k), 0};
    hsize_t count[2] = {1, 2};
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Create memory dataspace for just this row
    DataSpace memspace(2, count);

    // Write the row
    dataset.write(data, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
}



void Output::write_alpha_vs_time(int ts)
{
    std::string datasetName = "electronegativity";
    
    hsize_t ny = int(domain.NUM_TS / domain.write_interval) + 1; // rows (time steps)
    hsize_t dims[2] = {ny, 2}; // 2 columns: [time, alpha]

    // Create dataspace for full dataset
    DataSpace dataspace(2, dims);

    // Create or open the dataset
    H5::DataSet dataset;
    if (!time_group.exists(datasetName))
    {
        dataset = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    }
    else
    {
        dataset = time_group.openDataSet(datasetName);
    }

    // Prepare the data for this timestep
    int k = ts / domain.write_interval;
    double data[1][2];
    data[0][0] = ts * domain.DT;                    // Time
    data[0][1] = domain.electronegativity;   // Electronegativity value

    // Define hyperslab in the dataset where this row will be written
    hsize_t offset[2] = {static_cast<hsize_t>(k), 0};
    hsize_t count[2] = {1, 2};
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Create memory dataspace for just this row
    DataSpace memspace(2, count);

    // Write the row
    dataset.write(data, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
}



void Output::diagnostics(int ts, std::vector<Species> &species_list)
{
    
    if (domain.diagtype == "off")
    {
        std::cout << "TS: " << ts << std::endl;
        return;
    }

    auto norm_species = (domain.normscheme == 2 || domain.normscheme == 4) ? species_list[1] : species_list[0];
    
    
    if(domain.diagtype == "basic")
    {
        double max_phi = domain.phi(0);
        for(int i = 0; i < domain.ni; i++)
        {
            if (domain.phi(i) > max_phi)
            {
                max_phi = domain.phi(i);
            }
        }

        if(domain.SolverType == "direct" || domain.SolverType == "spectral")
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

        std::cout<<endl;
    }

    if(domain.diagtype == "full")
    {
        double max_phi = domain.phi(0);
        for(int i = 0; i < domain.ni; i++)
        {
            if (domain.phi(i) > max_phi)
            {
                max_phi = domain.phi(i);
            }
        }

        if(domain.SolverType == "direct" || domain.SolverType == "spectral")
        {
            std::cout << "TS: " << ts << " \t max_phi: " << std::fixed << std::setprecision(2) << (max_phi - domain.phi(0));
        }
        else
        {
            std::cout << "TS: " << ts << "\t" << "norm:" << domain.norm << " delta_phi: " << std::fixed << std::setprecision(2) << (max_phi - domain.phi(0));
        }

        if(domain.bc =="open")
        {
            //std::cout<<"\t"<<std::fixed << std::setprecision(precision)<<"rhoL|currentL:"<<domain.vL<<"|"<<domain.I_leftwall<<"\t"<<"rhoR|currentR:"<<domain.vR<<"|"<<domain.I_rightwall;
        }

        if(domain.bc == "open" || domain.enable_ionization_collision  == true || domain.enable_e_detach_collision == true)
        {
            for (Species& sp : species_list)
            {
                std::cout << " n_" << std::setw(4) << sp.name << ":" << sp.part_list.size();
            }
        }
        
        if(domain.enable_e_detach_collision == true)
        {
            std::cout << " alpha:" << domain.electronegativity;
        }

        if(domain.normscheme == 5)
        {
            std::cout << std::scientific << std::setprecision(precision);
        }
        else
        {
            std::cout << std::fixed << std::setprecision(precision);
        }

        double total_kinetic_energy = 0.0;
        vec<double> total_ke_components(3);
        double total_px = 0.0;
        double total_py = 0.0;
        double total_pz = 0.0;
        double potential_energy_val;

        if ((domain.bc == "pbc" || domain.bc == "rbc") && Energy_plot == 1)
        {
            potential_energy_val = domain.ComputePE(norm_species);
            for (Species& sp : species_list)
            {
                vec<double> ke = sp.Compute_KE(norm_species);
                total_ke_components += ke;
                total_kinetic_energy += ke(0) + ke(1) + ke(2);
                
                //std::cout << " KE_" << std::setw(4) << sp.name << ": " << ke(0) + ke(1) + ke(2);
            }

            std::cout << " KE_: " << total_kinetic_energy;
            std::cout << " PE_: " << potential_energy_val;
            std::cout << " TE_: " << total_kinetic_energy + potential_energy_val;

            /*
            for (Species& sp : species_list)
            {
                //total_momentum += sp.Compute_Momentum(norm_species);
                auto [px, py, pz] = sp.Compute_Momentum(norm_species);
                total_px += px;
                total_py += py;
                total_pz += pz;
            }*/
            //std::cout << " px: " << total_px << " py: " << total_py << " pz: " << total_pz<<" P:"<<sqrt(total_px*total_px + total_py*total_py + total_pz*total_pz);
            std::cout<< " delta_g:"<<domain.delta_g;
            double average_coll_freq = domain.avg_coll_freq/domain.W;
            //std::cout<< " \u03BD: "<< average_coll_freq;
            nu_avg.push_back(average_coll_freq);

            time_steps.push_back(static_cast<double>(ts));
            kinetic_energy.push_back(total_kinetic_energy);
            potential_energy.push_back(potential_energy_val);
            total_energy.push_back(total_kinetic_energy + potential_energy_val);
            Ke_x.push_back(total_ke_components(0));
            Ke_y.push_back(total_ke_components(1));
            Ke_z.push_back(total_ke_components(2));
            
        }
        std::cout << std::endl;

        //Plotting with matplotlibcpp
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

        //time.push_back(ts * domain.DT);

        if (Energy_plot == 1)
        {
            ke.push_back(total_kinetic_energy);
            pe.push_back(potential_energy_val);
            total_energie.push_back(total_kinetic_energy + potential_energy_val);
        }

        //if (domain.bc == "open")
        //{
        //    num1.push_back(species_list[0].part_list.size());
        //    num2.push_back(species_list[1].part_list.size());
        //}

        if (phase_plot == 1 )
        {
            //x.clear();
            //vx.clear();
            for(auto &part : species_list[species_index].part_list)
            {
                x.push_back(part.x);
                vx.push_back(part.vx);
            }
        }

        if (Chargedensity_plot == 1)
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

        if((Potentialfield_plot == 1 || domain.bc == "open") )
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

        if(dft_flag == 1 && domain.SolverType == "spectral")
        {
            k.clear();
            rho_k.clear();
            for(int i = 0; i < (domain.ni/2) +1 ; i++)
            {
                k.push_back(domain.dft_k(i));
                rho_k.push_back(domain.dft_value(i));
            }
        }

        plt::ion(); 

        if (Energy_plot == 1 && (domain.bc == "pbc" || domain.bc == "rbc"))
        {
            plt::figure(1);
            plt::clf();
            if (keflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full")
            {
                plt::named_plot("kinetic energy", time_steps, kinetic_energy , "r-");
            }
            if (peflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc") && domain.diagtype == "full")
            {
                plt::named_plot("potential energy", time_steps, potential_energy, "b-");
            }
            if (teflag == 1 && (domain.bc == "pbc" || domain.bc == "rbc" ) && domain.diagtype == "full")
            {
                plt::named_plot("total energy", time_steps, total_energy, "g-");
            }
            plt::xlabel("Time");
            plt::ylabel("Energy");
            plt::legend();
        }


        double marker_size = 1.0;
        std::map<std::string, std::string> style_options1 = {{"color", "blue"}, {"marker", "o"}};   
        std::map<std::string, std::string> style_options2 = {{"color", "red"}, {"marker", "o"}};

        if (phase_plot == 1)
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
            plt::xlim(0,domain.ni);
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


        if (Chargedensity_plot == 1)
        {   
            plt::figure(3);
            plt::clf();
            plt::named_plot("charge-density", lenght, charge_density, "r-");
            plt::xlabel("x");
            plt::ylabel("rho");
            plt::legend(); 
        }

        if((Potentialfield_plot == 1 || domain.bc == "open"))
        {
            plt::figure(4);
            plt::clf();
            plt::named_plot("potential", lenght, pot, "r-");
            //plt::named_plot("Electricfield", lenght, efield, "b-");
            plt::xlabel("x");
            //plt::ylim(-1000,1000);
            plt::ylabel("phi/Efield");
            plt::legend(); 
        }

        /*
        if(domain.bc == "open")
        {   
            plt::figure(5);
            plt::clf();
            plt::named_plot("electron", time_steps, num1, "r-");
            plt::named_plot("ion", time_steps, num2, "b-");
            //plt::named_plot("Electricfield", lenght, efield, "b-");
            plt::xlabel("time");
            plt::ylabel("simulation particle");
            plt::legend();

        }*/

       
        if(dft_flag == 1 && domain.SolverType == "spectral")
        {
            plt::figure(6);
            plt::clf();
            plt::named_plot("fourier transformed charge density", k, rho_k, "r-");
            plt::xlabel("k");
            plt::ylabel("rho_k");
            plt::legend(); 
        }


        if (ke_components == 1)
        {
            plt::figure(7);
            plt::clf();
            plt::plot(time_steps, Ke_x, {{"label", "KE_x"}, {"color", "blue"}});
            plt::plot(time_steps, Ke_y, {{"label", "KE_y"}, {"color", "red"}});
            plt::plot(time_steps, Ke_z, {{"label", "KE_z"}, {"color", "green"}});
            plt::title("Kinetic Energy Components vs Time");
            plt::xlabel("Time Step");
            plt::ylabel("Kinetic Energy component");
            plt::legend({{"loc", "upper right"}});
        }
        
        if(coll_freq_plot == 1)
        {
            plt::figure(8);
            plt::clf();
            plt::named_plot("average collision frequency", time_steps, nu_avg, "r-");
            plt::title("Average Collision Frequency vs Time");
            plt::xlabel("Time Step");
            plt::ylabel("\u03BD");
            plt::legend({{"loc", "upper right"}});
        }
        
        plt::pause(0.1);

        // Show plot
        plt::show();
    }
}



void Output::write_extra(double z1, double z2, double z3)
{
    // Create a subgroup under metadata or fielddata
    Group extra_group = file.createGroup("/extra_results");

    // Write attributes
    extra_group.createAttribute("z1", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &z1);

    extra_group.createAttribute("z2", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &z2);

    extra_group.createAttribute("z3", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &z3);

    extra_group.close();
}
