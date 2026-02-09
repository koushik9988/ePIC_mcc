#include "species.h"


/* note :
"-> " symbol means dereferencing a pointer , here we created pointers to instances of class named domain.
so to dereference the instances we have to use the symbol "->"*/

Species::Species(string name, double defaultmass, double charge, double defaultspwt, double temp, int numparticle, double vs, double fract_den , std:: string initialization_pos, string initialization_vel, Domain &domain):domain(domain)
{
    this->name = name;
    this->defaultmass = defaultmass;
    this->charge = charge;
    this->defaultspwt = defaultspwt;
    this->temp = temp;
    this->numparticle = numparticle;
    this->vs = vs;
    this->fract_den = fract_den;
    this->initialization_pos = initialization_pos;
    this->initialization_vel = initialization_vel;

    den = vec<double>(domain.ni);
    coll_rate = vec<double>(domain.ni);

    velmesh = vec<double>(domain.ni);
    power_rate = vec<double>(domain.ni);
    current_density = vec<double>(domain.ni);  
    
    //velprev = vec<double>(numparticle);

    //for(int i = 0; i < domain.num_threads; i++)
    //{
    //    domain.buffers.emplace_back(domain.ni);
    //}
    
}

void Species::AddParticle(Particle part)
{
    if (std::isnan(part.mass))
    {
        part.mass = this->defaultmass;
    }
    if (std::isnan(part.spwt))
    {
        part.spwt = this->defaultspwt;
    }

    this->part_list.push_back(part);
}


void Species::Push_species_serial(vector<Species> &species_list, int sub_cycle)
{
    double dt = domain.DT * sub_cycle;

    for (Particle &part : part_list)
    {
        double lc = part.x;
        double qm = charge / part.mass;

        double part_ef = domain.Interpolate(lc, domain.ef);

        // half step E
        double accel = qm * ((domain.density*Const::QE*domain.L) / (Const::EPS_0*domain.W*domain.vel_norm)) * part_ef;

        double vx_minus = part.vx + 0.5 * accel * dt;
        double vy_minus = part.vy;
        double vz_minus = part.vz;

        double vx_plus, vy_plus, vz_plus;

        // Boris push magnetic field       
        if (domain.B != 0.0)
        {
            vx_minus *= domain.vel_norm;
            vy_minus *= domain.vel_norm;
            vz_minus *= domain.vel_norm;

            double Bx = domain.B * sin(domain.theta) * cos(domain.azimuth);
            double By = domain.B * sin(domain.theta) * sin(domain.azimuth);
            double Bz = domain.B * cos(domain.theta);

            // t = (q/m) * B * dt / 2
            double tx = 0.5 * qm * Bx * (dt / domain.W);
            double ty = 0.5 * qm * By * (dt / domain.W);
            double tz = 0.5 * qm * Bz * (dt / domain.W);

            // v_prime = v_minus + v_minus x t
            double vpx = vx_minus + (vy_minus*tz - vz_minus*ty);
            double vpy = vy_minus + (vz_minus*tx - vx_minus*tz);
            double vpz = vz_minus + (vx_minus*ty - vy_minus*tx);

            // s = 2 / (1 + |t|^2)
            double s = 2.0 / (1.0 + tx*tx + ty*ty + tz*tz);

            // v_plus = v_minus + v_prime x (s * t)
            // sx = s * tx  | sy = s * ty |  sz = s * tz
            vx_plus = vx_minus + (vpy*(s*tz) - vpz*(s*ty));
            vy_plus = vy_minus + (vpz*(s*tx) - vpx*(s*tz));
            vz_plus = vz_minus + (vpx*(s*ty) - vpy*(s*tx));

            vx_plus /= domain.vel_norm;
            vy_plus /= domain.vel_norm;
            vz_plus /= domain.vel_norm;
        }
        else
        {
            vx_plus = vx_minus;
            vy_plus = vy_minus;
            vz_plus = vz_minus;
        }

        //Half-step acceleration by E field 
        part.vx = vx_plus + 0.5 * accel * dt;
        part.vy = vy_plus;
        part.vz = vz_plus;

        
        double old_pos = part.x;
        part.x += ((domain.vel_norm)/(domain.L*domain.W)) * part.vx * dt;

        double new_pos = part.x;

        
        if (fabs(new_pos - old_pos) >= 1)
        {
            if (name == "electron")
            {
                domain.ele_cross++;
                domain.crossed_cellno_ele += int(fabs(new_pos - old_pos));
            }
            else if (name == "ion")
            {
                domain.ion_cross++;
                domain.crossed_cellno_ion += int(fabs(new_pos - old_pos));
            }
        }

        // boundary conditions
        if (domain.bc == "pbc")
        {
            if (part.x < domain.x0)
            {
                part.x += domain.xL;
            }
            else if (part.x >= domain.x0 + domain.xL)
            {
                part.x -= domain.xL;
            }
        }
        else if (domain.bc == "rbc")
        {
            if (part.x < domain.x0)
            {
                part.x = domain.x0 + (domain.x0 - part.x);
                part.vx = -part.vx;
            }
            else if (part.x >= domain.x0 + domain.xL)
            {
                part.x = domain.x0 + domain.xL - (part.x - (domain.x0 + domain.xL));
                part.vx = -part.vx;
            }
        }
        else if (domain.bc == "open")
        {
            if (part.x < domain.x0 || part.x >= domain.x0 + domain.xL)
            {
                part.alive = false; 

                if (part.x < domain.x0)
                {
                    domain.wall_left = true;
                    if (IsIon()) SEE(species_list[0], domain);
                    domain.vL += (IsIon()?1:-1) * part.spwt/domain.density;
                }
                else
                {
                    domain.wall_left = false;
                    if (IsIon()) SEE(species_list[0], domain);
                    domain.vR += (IsIon()?1:-1) * part.spwt/domain.density;
                }
            }
        }
    }

    if (domain.bc == "open")
    {
        
        /*
            remove_if function doesnot actually erase elements from the vector, it just moves the "alive" particles to the front and returns an iterator to the new end of the vector. 
            Then the erase function is called to remove the "dead" particles from the vector based on the new end iterator returned by remove_if. This approach is efficient and
            avoids issues with modifying the vector while iterating over it in a multithreaded context.
            Argument of remove_if is iterator to the beginning of the vector, iterator to the end of the vector and a condition defined by a
            lambda function that returns true for particles that are not alive (i.e., should be removed).
        */
        
        part_list.erase(std::remove_if(part_list.begin(), part_list.end(),[](const Particle& p){ return !p.alive; }),part_list.end());

        domain.I_leftwall  = domain.vL / dt;
        domain.I_rightwall = domain.vR / dt;
    }
}

void Species::update(int start, int end, int sub_cycle, vector<Species> &species_list)
{
    //Use iterators for parallel processing
    auto it = part_list.begin() + start;
    auto end_it = part_list.begin() + end;
    
    int index = 0;

    while (it != end_it)
    {
        Particle &part = *it;
        // double lc = domain.XtoL(part.pos);
        double lc = part.x;
        // Compute charge to mass ratio
        double qm = charge / part.mass;

        // Interpolate electric field at particle location
        double part_ef = domain.Interpolate(lc, domain.ef);

        double dt = domain.DT * sub_cycle;
        
        // --- Half-step acceleration by E field
        double vx_minus = part.vx  +  0.5* qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_ef*(domain.DT*sub_cycle);
        double vy_minus = part.vy;
        double vz_minus = part.vz;

        double vx_plus, vy_plus, vz_plus;

        if (domain.B != 0.0 )
        {
            // Unnormalize velocities
            vx_minus *= domain.vel_norm;
            vy_minus *= domain.vel_norm;
            vz_minus *= domain.vel_norm;

            double B = domain.B;
            double theta = domain.theta;
            double azimuth = domain.azimuth;

            double Bx = B * sin(theta) * cos(azimuth);
            double By = B * sin(theta) * sin(azimuth); 
            double Bz = B * cos(theta);

            // t = (q/m) * B * dt / 2
            double tx = 0.5 * qm * Bx * (dt / domain.W);
            double ty = 0.5 * qm * By * (dt / domain.W);
            double tz = 0.5 * qm * Bz * (dt / domain.W);

            // v_prime = v_minus + v_minus x t
            double v_prime_x = vx_minus + (vy_minus * tz - vz_minus * ty);
            double v_prime_y = vy_minus + (vz_minus * tx - vx_minus * tz);
            double v_prime_z = vz_minus + (vx_minus * ty - vy_minus * tx);

            // s = 2 / (1 + |t|^2)
            double t2 = tx * tx + ty * ty + tz * tz;
            double s = 2.0 / (1.0 + t2);

            // v_plus = v_minus + v_prime x (s * t)
            // sx = s * tx  | sy = s * ty |  sz = s * tz
            vx_plus = vx_minus + (v_prime_y * (s * tz) - v_prime_z * (s * ty));
            vy_plus = vy_minus + (v_prime_z * (s * tx) - v_prime_x * (s * tz));
            vz_plus = vz_minus + (v_prime_x * (s * ty) - v_prime_y * (s * tx));

            // Normalize velocities
            vx_plus /= domain.vel_norm;
            vy_plus /= domain.vel_norm;
            vz_plus /= domain.vel_norm;
        }
        else
        {
            vx_plus = vx_minus;
            vy_plus = vy_minus;
            vz_plus = vz_minus;
        }

        // --- Half-step acceleration by E field (again)
        part.vx = vx_plus + 0.5* qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_ef*(domain.DT*sub_cycle);
        part.vy = vy_plus;
        part.vz = vz_plus;

        double old_pos = part.x;
        // --- Position update
        part.x += ((domain.vel_norm)/(domain.L*domain.W))*part.vx*(domain.DT*sub_cycle);

        double new_pos = part.x;
       
        // Check for particle crossing
        if (fabs(new_pos - old_pos) >= 1)
        {
            if (name == "electron")
            {
                domain.ele_cross++;
                domain.crossed_cellno_ele += int(fabs(new_pos - old_pos));
            }
            if (name == "ion")
            {
                domain.ion_cross++;
                domain.crossed_cellno_ion += int(fabs(new_pos - old_pos));
            }
            // cerr << "Warning: Particle crossing a full cell!" << endl;
        }

        // Handle boundary conditions
        if (domain.bc == "pbc")
        {
            if (part.x < domain.x0)
            {
                part.x += domain.xL;
            }
            else if (part.x >= domain.x0 + domain.xL)
            {
                part.x -= domain.xL;
            }
        }

        else if (domain.bc == "rbc")
        {
            if (part.x < domain.x0)
            {
                part.x = domain.x0 + (domain.x0 - part.x); // Reflect position
                part.vx = -part.vx; // Reverse velocity
            }
            else if (part.x >= domain.x0 + domain.xL)
            {
                part.x = domain.x0 + domain.xL - (part.x - (domain.x0 + domain.xL)); // Reflect position
                part.vx = -part.vx; // Reverse velocity
            }
        }

        else if(domain.bc == "open")
        {
            if (part.x < domain.x0 || part.x >= domain.x0 + domain.xL)
            {
                // Flag particle as dead instead of erasing (thread-safe)
                part.alive = false;
                
                if(part.x < domain.x0)
                {
                    domain.wall_left = true;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }

                    // Accumulate charge at the left wall
                    if (IsIon())
                    {
                        domain.vL += 1*part.spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vL += -1*part.spwt/domain.density ; 
                    }
                }
                if(part.x >= domain.x0 + domain.xL)
                {
                    domain.wall_left = false;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }
                    //new code line
                    // Accumulate charge at the right wall
                    if (IsIon())
                    {
                        domain.vR += 1*part.spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vR += -1*part.spwt/domain.density; 
                    }
                }
            }
        }
        ++it;
    }
}

void Species::Push_species(vector<Species> &species_list, int sub_cycle)
{
    
    if(domain.push_parallal)
    {
        int particles_per_threads = numparticle/domain.num_threads;
       
        std::vector<thread> threads;

        for (int i = 0; i < domain.num_threads; ++i) 
        {
            int start = i * particles_per_threads;
            int end = (i == domain.num_threads - 1) ? part_list.size() : (i + 1) * particles_per_threads;

            threads.emplace_back(&Species::update, this, start, end, sub_cycle, std::ref(species_list));
        }

        // Join threads to wait for their completion
        for (auto &thread : threads) 
        {
            thread.join();
        }
        
        // After all threads complete, remove dead particles (thread-safe)
        if(domain.bc == "open")
        {
            /*
            remove_if function doesnot actually erase elements from the vector, it just moves the "alive" particles to the front and returns an iterator to the new end of the vector. 
            Then the erase function is called to remove the "dead" particles from the vector based on the new end iterator returned by remove_if. This approach is efficient and
            avoids issues with modifying the vector while iterating over it in a multithreaded context.
            Argument of remove_if is iterator to the beginning of the vector, iterator to the end of the vector and a condition defined by a
            lambda function that returns true for particles that are not alive (i.e., should be removed).
            */
        
            part_list.erase(std::remove_if(part_list.begin(), part_list.end(),[](const Particle& p) { return !p.alive; }),part_list.end());
              
            // Compute wall currents after charge accumulation
            double J_L = domain.vL / (domain.DT * sub_cycle);
            double J_R = domain.vR / (domain.DT * sub_cycle);
            
            domain.I_leftwall = J_L;
            domain.I_rightwall = J_R;
        }
    }
    else
    {
        Push_species_serial( species_list, sub_cycle);
    }
}

void Species::ScatterSpecies_serial()
{   
    //reset density
    den = 0;
    for(Particle &part : part_list)
    {
        double lc = domain.XtoL(part.x);
        //double lc = part.pos;
        domain.Scatter(lc,part.spwt,den);
    }

    den /= (domain.dx*domain.L);
    den /= domain.density;

    if(domain.bc == "pbc")
    {
        den(0) += den(domain.ni-1);
        den(domain.ni-1) =  den(0);
    }
    else if(domain.bc == "open" || domain.bc == "rbc")
    {
        den(0) *= 2;
        den(domain.ni-1) *= 2;
    }
}

/// mesh averaged velocity
void Species::ScatterVel_serial()
{   
    //reset density
    velmesh = 0;
    for(Particle &part : part_list)
    {
        double lc = domain.XtoL(part.x);
        domain.Scatter(lc,part.spwt*part.vx,velmesh);
    }

    velmesh /= (domain.dx*domain.L);

    for(int i = 0 ; i < domain.ni ; i++)
    {
        velmesh(i) /= den(i)*domain.density;
    }
    //velmesh /= den*domain.density;

    if(domain.bc == "pbc")
    {
        velmesh(0) += velmesh(domain.ni-1);
        velmesh(domain.ni-1) =  velmesh(0);
    }
    else if(domain.bc == "open" || domain.bc == "rbc")
    {
        velmesh(0) *= 2;
        velmesh(domain.ni-1) *= 2;
    }
}

void Species::parall_deposit(int threadidx, int start, int end)
{
    //domain.buffers[threadidx].clear();
    for(int i = start; i < end; i++)
    {
        Particle &part = part_list[i];
        double lc = domain.XtoL(part.x);
        domain.Scatter(lc,part.spwt,local_buffers[threadidx]);
    }
}

void Species::ScatterSpecies()
{
    if(domain.deposit_parallal)
    {
        int particles_per_threads = part_list.size()/domain.num_threads;

        ////
        local_buffers.clear();
        local_buffers.resize(domain.num_threads);
        
        // Explicitly zero all elements in each buffer to avoid uninitialized memory issues (Is binary addition associative and commutative ?)
        for (int i = 0; i < domain.num_threads; i++) 
        {
            local_buffers[i] = vec<double>(domain.ni);
            for (int j = 0; j < domain.ni; j++)
            {
                local_buffers[i](j) = 0.0;
            }
        }
        ////       
        std::vector<thread> threads;
        
        for (int i = 0; i < domain.num_threads; i++) 
        {
            int start = i * particles_per_threads;
            int end = (i == domain.num_threads - 1) ? part_list.size() : start + particles_per_threads;
            threads.emplace_back(&Species::parall_deposit,this,i, start, end);
        }

        //Join threads 
        for (auto &thread : threads) 
        {
            thread.join();
        }
        
        //reset density
        den.clear();

        for(int i = 0 ; i < domain.num_threads; i++)
        {
            den += local_buffers[i];
        }

        den /= (domain.dx*domain.L);
        den /= domain.density;

        if(domain.bc == "pbc")
        {
            den(0) += den(domain.ni-1);
            den(domain.ni-1) =  den(0);
        }
        else if(domain.bc == "open" || domain.bc == "rbc")
        {
            den(0) *= 2;
            den(domain.ni-1) *= 2;
        }
    }
    else
    {
        ScatterSpecies_serial();
    }
}


void Species::Rewind_species()
{
    for (Particle &part: part_list)
    {
        double lc = domain.XtoL(part.x);
        //double lc = part.pos;
        // compute charge to mass ratio
        double qm = charge/part.mass;
        //cout<<qm<<endl;
        double part_ef = domain.Interpolate(lc,domain.ef);
        double wl = domain.LDe*domain.LDe*domain.wpe*domain.wpe;

        part.vx -= 0.5*qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_ef*domain.DT;
                
    }
}

vec<double> Species::Compute_KE(Species &species)
{
    vec<double> ke(3);
    for (Particle &part:part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke(0) += (part.vx*part.vx)*(domain.vel_norm)*(domain.vel_norm)*part.mass*part.spwt;
        ke(1) += (part.vy*part.vy)*(domain.vel_norm)*(domain.vel_norm)*part.mass*part.spwt;
        ke(2) += (part.vz*part.vz)*(domain.vel_norm)*(domain.vel_norm)*part.mass*part.spwt;
    }
    /*Multiply 0.5*mass for all particles*/
    ke(0) *= 0.5;
    ke(1) *= 0.5;
    ke(2) *= 0.5;
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (species.temp*Const::eV)*(species.defaultspwt)*species.numparticle;
    if(domain.normscheme == 5)
    {Th = 1;}

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke(0) = ke(0)/Th;
    ke(1) = ke(1)/Th;
    ke(2) = ke(2)/Th;
    return ke;
}

std::tuple<double, double, double> Species::Compute_Momentum(Species &species)
{
    double p_x = 0.0, p_y = 0.0, p_z = 0.0;

    for (Particle &part : part_list)
    {
        // Compute total momentum in each direction (average velocity)
        p_x += (part.vx) * domain.vel_norm *part.mass*part.spwt;
        p_y += (part.vy)* domain.vel_norm *part.mass*part.spwt;
        p_z += (part.vz)* domain.vel_norm *part.mass*part.spwt;
    }

    
    // Normalize momentum
    double Thp;
    if (domain.normscheme == 5)
    {
        Thp = sqrt(domain.energy_scale * species.defaultspwt * species.numparticle * Const::ME);
    }
    else
    {
        Thp = sqrt((species.temp * Const::eV) / Const::ME) * species.defaultspwt * species.numparticle * Const::ME;
    }

    p_x /= Thp;
    p_y /= Thp;
    p_z /= Thp;

    return std::make_tuple(p_x, p_y, p_z);
}

void SEE(Species &species, Domain &domain)
{
    int num_see = int(domain.see_rate + rnd());
    if(domain.see_rate == 0)
    {
        num_see = 0;
    }
    double x,v;
    for(int i = 0; i < num_see; i++)
    {
        if(domain.wall_left)
        {
            x = domain.x0;
            v = Init::SampleVel(species,domain.tempwall);
        }
        else
        {
            x = domain.xL-rnd();
            v = - Init::SampleVel(species,domain.tempwall);
        }

        species.AddParticle(Particle(x,v,0,0,1));
    }
}

bool Species::IsIon()
{
    if(charge > 0)
    {return true;}
    else
    {return false;}
}


void Species::ComputeCurrentDensity()
{
    //Reset current density to zero
    current_density = 0;
    
    //Calculate current density at each grid point: j = q * n * v
    //where:
    //q = charge of the species
    //n = number density at grid point (den)
    //v = mean velocity at grid point (velmesh)
    
    for (int i = 0; i < domain.ni; i++)
    {
        // Current density = charge * density * velocity * unnorm 
        current_density(i) = charge * den(i) * velmesh(i) * domain.density * domain.vel_norm;
    }
}

void Species::ComputeHeatingRate()
{
    // Reset heating rate to zero
    power_rate = 0;
    
    //Calculate heating rate at each grid point: P = j * E
    //where:
    //j = current density at grid point (current_density)
    //E = electric field at grid point (domain.ef)
    
    for (int i = 0; i < domain.ni; i++)
    {
        // Heating rate (power absorption) = current density * electric field
        power_rate(i) = current_density(i) * domain.ef(i) * ((domain.density * Const::QE * domain.L) / Const::EPS_0);
    }
}