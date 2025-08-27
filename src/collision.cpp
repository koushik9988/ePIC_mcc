/*
Implementation file for MCC collisional models
*/

#include "collision.h"

CollisionHandler::CollisionHandler(Domain &domain, int gas, const char* data_path) : domain(domain), gas(gas)
{
    sigma_tot_e.resize(CS_RANGES, 0.0);
    sigma_tot_pion.resize(CS_RANGES, 0.0);
    sigma_tot_nion.resize(CS_RANGES, 0.0);

    sigma_ela.resize(CS_RANGES, 0.0);
    sigma_exc.resize(CS_RANGES, 0.0);
    sigma_ionz.resize(CS_RANGES, 0.0);
    sigma_mex_pion.resize(CS_RANGES, 0.0);
    sigma_e_detach.resize(CS_RANGES,0.0); 

    if (gas == HYDROGEN)
    {
        setup_hydrogen_cross_sections(data_path);
    }
    // Argon doesn't need setup (uses analytical expressions)
}

CollisionHandler::CollisionHandler(Domain &domain, int gas) : domain(domain), gas(gas)
{
    sigma_tot_e.resize(CS_RANGES, 0.0);
    sigma_tot_pion.resize(CS_RANGES, 0.0);
    sigma_tot_nion.resize(CS_RANGES, 0.0);

    sigma_ela.resize(CS_RANGES, 0.0);
    sigma_exc.resize(CS_RANGES, 0.0);
    sigma_ionz.resize(CS_RANGES, 0.0);
    sigma_mex_pion.resize(CS_RANGES, 0.0);
    sigma_e_detach.resize(CS_RANGES,0.0); 
}


void CollisionHandler::set_electron_cross_sections()
{
    std::vector<double> energy(CS_RANGES);
    energy[0] = DE_CS;
    std::generate(energy.begin() + 1, energy.end(), [i = 1]() mutable { return DE_CS * (i++); });

    // Select gas-specific cross-section functions
    if (gas == ARGON)
    {
        std::transform(energy.begin(), energy.end(), sigma_ela.begin(),[this](double energy_val) { return compute_elastic_CS_ar(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_exc.begin(),[this](double energy_val) { return compute_excitation_CS_ar(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_ionz.begin(),[this](double energy_val) { return compute_ionization_CS_ar(energy_val, domain); });
    }
    else if (gas == HYDROGEN)
    {
        //display::print("executed");
        std::transform(energy.begin(), energy.end(), sigma_ela.begin(),[this](double energy_val) { return compute_elastic_CS_h(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_exc.begin(),[this](double energy_val) { return compute_excitation_CS_h(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_ionz.begin(),[this](double energy_val) { return compute_ionization_CS_h(energy_val, domain); });
    }
    else
    {
        std::cerr << "Error: Unknown gas type\n";
    }
}

void CollisionHandler::set_pion_cross_sections()
{

    std::vector<double> energy(CS_RANGES);
    energy[0] = DE_CS;
    std::generate(energy.begin() + 1, energy.end(), [i = 1]() mutable { return DE_CS * (i++); });

    if (gas == HYDROGEN)
    {
        std::transform(energy.begin(), energy.end(), sigma_mex_pion.begin(),[this](double energy_val) { return compute_mex_cs_pion(energy_val, domain); });
        std::transform(energy.begin(), energy.end(), sigma_e_detach.begin(),[this](double energy_val) { return compute_e_detach(energy_val, domain); });
    }
    else
    {
        std::cerr << "Error: Unknown gas type\n";
    }


}

void CollisionHandler::calc_total_cross_sections()
{
    for (size_t i = 0; i < CS_RANGES; ++i)
    {
        sigma_tot_e[i] = (sigma_ela[i] + sigma_exc[i] + sigma_ionz[i]) * domain.GAS_DENSITY;
        sigma_tot_pion[i] = (sigma_mex_pion[i]) * domain.GAS_DENSITY;
        sigma_tot_nion[i] = sigma_e_detach[i] * domain.GAS_DENSITY;
    }
}

double CollisionHandler::max_electron_coll_freq()
{
    double e, v, nu, nu_max = 0.0;
    for (int i = 0; i < CS_RANGES; ++i)
    {
        e = i * DE_CS;
        v = sqrt(2.0 * e * Const::eV / Const::ME);
        nu = v * sigma_tot_e[i];
        if (nu > nu_max)
        {
            nu_max = nu;
        }
    }
    return nu_max;
}

std::pair<double, double> CollisionHandler:: electron_min_max_mean_free_path()
{
    double e, v, sigma, lambda;
    double lambda_min = 1e10;  // a very large starting value
    double lambda_max = 0.0;   // a very small starting value

    for (int i = 0; i < CS_RANGES; ++i)
    {
        e = i * DE_CS;
        sigma = sigma_tot_e[i];
        //v = sqrt(2.0 * e * Const::eV / Const::ME); // m/s
        lambda = 1.0 / (sigma);       // m

        if (lambda < lambda_min) lambda_min = lambda;
        if (lambda > lambda_max) lambda_max = lambda;
    }

    return {lambda_min, lambda_max};
}


void CollisionHandler::collision_electron(double xe, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2)
{

    const double F1 = Const::ME / (Const::ME + species2.defaultmass);
    const double F2 = species2.defaultmass / (Const::ME + species2.defaultmass);
    double  t0, t1, t2;//rnd;
    double  g, g2, gx, gy, gz, wx, wy, wz, theta, phi;
    double  chi, eta, chi2, eta2, sc, cc, se, ce, st, ct, sp, cp, energy, e_sc, e_ej;
    double g_before, g_after;

    // Calculate relative velocity before collision & velocity of the center of mass before collision
    gx = vxe;
    gy = vye;
    gz = vze;
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    g_before = g;
    wx = F1 * vxe;
    wy = F1 * vye;
    wz = F1 * vze;

    //display::print("gbefore:",g); //gx,gy,gz giving nan value have to fix.

    // Find Euler angles
    if (gx == 0) 
    {
        theta = 0.5 * Const::PI;
    }
    else 
    {
        theta = atan2(sqrt(gy * gy + gz * gz),gx);
    }
    if (gy == 0)
    {
        if (gz > 0)
        {
            phi = 0.5 * Const::PI;
        }
        else
        {
            phi = - 0.5 * Const::PI;
        }
    }
    else
    {
        phi = atan2(gz, gy);
    }
    st  = sin(theta);
    ct  = cos(theta);
    sp  = sin(phi);
    cp  = cos(phi);

    // Choose the type of collision based on cross-sections
    t0 = sigma_ela[eindex];
    t1 = t0 + sigma_exc[eindex];
    t2 = t1 + sigma_ionz[eindex];
    //rnd = rand();  // Random number between 0 and 1

    if (rnd() < (t0 / t2))
    {  
        //display::print("elatic!");
        // Elastic scattering
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2*Const::PI * rnd();
        
    }
    else if (rnd() < (t1 / t2))
    {  
        //display::print("inelastic!");
        // Excitation
        energy = 0.5 * Const::ME *g * g * domain.vel_norm * domain.vel_norm;  // Energy in joules
        energy = fabs(energy - E_EXC_TH * Const::eV);  // Energy loss for excitation (exitation enegry converted to joule)
        g = sqrt(2.0 * energy / Const::ME);  // Relative velocity after energy loss
        g = g / domain.vel_norm;  // Normalize the velocity
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2*Const::PI * rnd();
    }
    else
    {  
        // Ionization
        //display::print("ionization!");
        energy = 0.5 * Const::ME * g * g * domain.vel_norm * domain.vel_norm; // Energy in joules
        energy = fabs(energy - E_ION_TH * Const::eV);  // Energy loss for ionization
        e_ej = 10.0 * tan(rnd() * atan(energy / Const::eV / 20.0)) * Const::eV;  // Emitted electron energy
        e_sc = fabs(energy - e_ej);  // Incoming electron energy after collision
        g = sqrt(2.0 * e_sc / Const::ME);
        g2 = sqrt(2.0 * e_ej / Const::ME);
        g = g / domain.vel_norm;  // Normalize the velocity
        g2 = g2 / domain.vel_norm;  // Normalize the velocity
        chi = acos(sqrt(e_sc / energy));  // Scattering angle for incoming electron
        chi2 = acos(sqrt(e_ej / energy));  // Scattering angle for emitted electron
        eta = 2*Const::PI * rnd();
        eta2 = eta + Const::PI;

        // Compute velocity of emitted electron
        sc = sin(chi2);
        cc = cos(chi2);
        se = sin(eta2);
        ce = cos(eta2);
        gx = g2 * (ct * cc - st * sc * ce);
        gy = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
        
        // Add the new electron and ion (normalize this as we have to unnormalize for previous calculations) (velocities are converted to lab gram from com frame here)
        species1.AddParticle(Particle(xe, (wx + F2 * gx), (wy + F2 * gy), (wz + F2 * gz), 0));
        species2.AddParticle(Particle(xe, Init::SampleVel(species2,species2.temp),  Init::SampleVel(species2,species2.temp),  Init::SampleVel(species2,species2.temp),0));
    }

    // Scatter the primary electron
    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);

    // Compute new relative velocity 
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // Post-collision velocity of the colliding electron in lab frame
    vxe = wx + F2 * gx;
    vye = wy + F2 * gy;
    vze = wz + F2 * gz;

    g_after = sqrt(vxe * vxe + vye * vye + vze * vze);

    domain.delta_g = (fabs(g_after - g_before));  // Change in velocity after collision
    //printf("g_before: %f, g_after: %f, delta_g: %f\n", g_before*domain.vel_norm, g_after*domain.vel_norm, domain.delta_g);

}


void CollisionHandler::collision_ion(double xe, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2,
     int eindex, Species &species1, Species &species2)
{
    double   g,gx,gy,gz,wx,wy,wz;
    double   theta,phi,chi,eta,st,ct,sp,cp,sc,cc,se,ce,t;
    
    // calculate relative velocity before collision
    // random Maxwellian target atom already selected (vx_2,vy_2,vz_2 velocity components of target atom come with the call)
    gx = vx_1-vx_2;
    gy = vy_1-vy_2;
    gz = vz_1-vz_2;
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    //in case of ion-neutral collision we have to consider neutral motion
    wx = (vx_1 * species1.defaultmass + vx_2 * species2.defaultmass)/(species1.defaultmass + species2.defaultmass);
    wy = (vy_1 * species1.defaultmass + vy_2 * species2.defaultmass)/(species1.defaultmass + species2.defaultmass);
    wz = (vz_1 * species1.defaultmass + vz_2 * species2.defaultmass)/(species1.defaultmass + species2.defaultmass);

    // find Euler angles:

    if (gx == 0)
    {
        theta = 0.5 * Const::PI;
    }
    else
    {
        theta = atan2(sqrt(gy * gy + gz * gz),gx);
    }
    if (gy == 0)
    {
        if (gz > 0)
        {
            phi = 0.5 * Const::PI;
        }
        else
        {
            phi = - 0.5 * Const::PI;
        }
    }
    else 
    {
        phi = atan2(gz, gy);
    }

    // determine the type of collision based on cross sections and generate scattering angle

    t  = sigma_mex_pion[eindex];
    
    if  (rnd() < 1) // 1 because we only considered only one type of process(elastic)
    {                                 
        chi = acos(1.0 - 2.0 * rnd()); 
        eta = 2 * Const::PI * rnd();                
    } 
                               
    sc  = sin(chi);
    cc  = cos(chi);
    se  = sin(eta);
    ce  = cos(eta);
    st  = sin(theta);
    ct  = cos(theta);
    sp  = sin(phi);
    cp  = cos(phi);

    // compute new relative velocity:

    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    //post-collision velocity of the ion

    double F1 = species2.defaultmass / (species1.defaultmass + species2.defaultmass);;
    vx_1 = wx + F1 * gx;
    vy_1 = wy + F1 * gy;
    vz_1 = wz + F1 * gz; 
}


void CollisionHandler::collision_detachment(double xe, double &vx_1, double &vy_1, double &vz_1, double &vx_2, double &vy_2, double &vz_2, 
    int eindex, Species &species1, Species &species2, Species &electron_species)
{
    // Relative velocity (g) and center-of-mass velocity (w)
    double gx = vx_1 - vx_2;
    double gy = vy_1 - vy_2;
    double gz = vz_1 - vz_2;
    double g_mag = sqrt(gx * gx + gy * gy + gz * gz);

    double m1 = species1.defaultmass;
    double m2 = species2.defaultmass;
    double m_total = m1 + m2;
    
    double wx = (vx_1 * m1 + vx_2 * m2) / m_total;
    double wy = (vy_1 * m1 + vy_2 * m2) / m_total;
    double wz = (vz_1 * m1 + vz_2 * m2) / m_total;

    double E_th_det = 2.25 * Const::eV; // Detachment threshold energy in Joules
    
    double E_kin_com = 0.5 * ((m1*m2)/m_total) * (g_mag * g_mag * domain.vel_norm * domain.vel_norm); //reduced mass in COM frame

    if (E_kin_com < E_th_det) return; // Not enough energy, no collision

    double E_available = E_kin_com - E_th_det;

    // Simple model for energy sharing
    double E_electron_ej = 0;//E_available * (1.0 - sqrt(rnd()));

    double v_e_com = sqrt(2.0 * E_electron_ej / Const::ME);

    // Isotropic scattering: generate random angles for the new particles.
    double chi = acos(1.0 - 2.0 * rnd()); 
    double eta = 2.0 * Const::PI * rnd();

    double sc = sin(chi);
    double cc = cos(chi);
    double se = sin(eta);
    double ce = cos(eta);
 
    // Velocity components for the new electron in COM frame
    double ve_com_x = v_e_com * sc * ce;
    double ve_com_y = v_e_com * sc * se;
    double ve_com_z = v_e_com * cc;

    electron_species.AddParticle(Particle(xe, (wx + ve_com_x / domain.vel_norm),(wy + ve_com_y / domain.vel_norm),(wz + ve_com_z / domain.vel_norm),1,NAN,species1.defaultspwt));
    //display::print("executed detachment collision");
    
}


//  (spwt?????)
/*
void CollisionHandler::handle_collisions(Species &projectile, Species &target_gas, Species &electron_species)
{
    for (Particle &part : projectile.part_list)
    {
        double v_sqr = (part.vx * part.vx + part.vy * part.vy + part.vz * part.vz) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);
        double energy = (0.5 * projectile.defaultmass * v_sqr) / Const::eV;
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        double sigma_tot = 0.0;

        if (projectile.name == "electron")
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion")
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }
        else if (projectile.name == "beam" )//&& projectile.charge == -1)
        {
            sigma_tot = sigma_tot_pion[energy_index];//sigma_tot_e[energy_index];
        }
        else
        {
            std::cerr << "CollisionHandler: Unknown projectile type: " << projectile.name << std::endl;
            continue;
        }

        double nu = sigma_tot * velocity;

        double p_coll = 0.5;//1 - exp(-nu * (domain.DT / domain.W));

        if (p_coll > rnd())
        {
            part.id  = 1; // Mark particle as colliding
            //domain.N_pcoll = domain.N_pcoll + 1;
            // Sample target gas particle velocity from a maxwellian distribution
            double vxg = Init::SampleVel(target_gas, target_gas.temp);
            double vyg = Init::SampleVel(target_gas, target_gas.temp);
            double vzg = Init::SampleVel(target_gas, target_gas.temp);

            if (projectile.name == "electron")
            {
                display::print("executed electron collision");
                collision_electron(part.x, part.vx, part.vy, part.vz,energy_index, projectile, target_gas);
                domain.N_pcoll++;
                //collision_electron(part.x, part.vx, part.vy, part.vz, energy_index, electron, target_gas);
            }
            else if (projectile.name == "ion")
            {
                display::print("executed ion collision");
                collision_ion(part.x, part.vx, part.vy, part.vz, vxg, vyg, vzg, energy_index, projectile, target_gas);
            }

            else if (projectile.name == "beam" && projectile.charge == -1)
            {
                display::print("executed beam collision");
                collision_detachment(part.x, part.vx, part.vy, part.vz,vxg, vyg, vzg,energy_index, projectile, target_gas, electron_species);
                part.id = -1; // Mark beam particle for deletion after collision
        
            }
        }
    }

    // Remove particles marked for deletion
    projectile.part_list.erase(std::remove_if(projectile.part_list.begin(), projectile.part_list.end(),[](const Particle &p) { return p.id == -1; }),projectile.part_list.end());
}*/



void CollisionHandler::handle_collisions(Species &projectile, Species &target_gas, Species &electron_species)
{
    for (auto it = projectile.part_list.begin(); it != projectile.part_list.end(); )
    {
        Particle &part = *it;

        double v_sqr = (part.vx * part.vx + part.vy * part.vy + part.vz * part.vz) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);
        double energy   = (0.5 * projectile.defaultmass * v_sqr) / Const::eV;
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        double sigma_tot = 0.0;
        if (projectile.name == "electron" || projectile.name == "ebeam") //name are hardcoded but may be i can use electron mass and charge as a condition
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion")
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }
        else if (projectile.name == "beam" || projectile.name == "negion") // && projectile.charge == -1
        {
            sigma_tot = sigma_tot_nion[energy_index]; // or sigma_tot_e if desired
        }
        

        double nu = sigma_tot * velocity;
        double p_coll = 1 - exp(-nu * (domain.DT / domain.W));

        if (p_coll > rnd())
        {
            it->id = 1;
            // Maxwellian sample for target velocity
            double vxg = Init::SampleVel(target_gas, target_gas.temp)/ domain.vel_norm; // normalize to lab frame
            double vyg = Init::SampleVel(target_gas, target_gas.temp)/ domain.vel_norm;; //normalize 
            double vzg = Init::SampleVel(target_gas, target_gas.temp) / domain.vel_norm;

            if (projectile.name == "electron" || projectile.name == "ebeam") // species name are hardcoded
            {
                //display::print("executed electron collision");
                collision_electron(part.x, part.vx, part.vy, part.vz, energy_index, projectile, target_gas);
                domain.N_ecoll++;
            }
            else if (projectile.name == "ion")
            {
                //display::print("executed ion collision");
                collision_ion(part.x, part.vx, part.vy, part.vz, vxg, vyg, vzg, energy_index, projectile, target_gas);
                domain.N_ioncoll++;
            }
            else if (projectile.name == "beam" || projectile.name == "negion")
            {
                //display::print("executed beam collision");
                collision_detachment(part.x, part.vx, part.vy, part.vz,vxg, vyg, vzg, energy_index, projectile, target_gas, electron_species);
                // erase this particle immediately
                it = projectile.part_list.erase(it);
                continue; // skip ++it, since erase returns next element
                domain.N_negioncoll++;
            }
        }

        ++it; // move to next particle if not erased
    }
}

double CollisionHandler::average_collision_frequency(Species &projectile)
{
    double total_nu = 0.0;
    int num_particles = projectile.part_list.size();

    for (const Particle &part : projectile.part_list)
    {
        // Calculate particle speed (velocity magnitude)
        double v_sqr = (part.vx * part.vx + part.vy * part.vy + part.vz * part.vz) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);

        // Energy in eV using mass of projectile species
        double energy = (0.5 * projectile.defaultmass * v_sqr) / Const::eV;

        // Get energy index
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);
        //display::print(energy_index);

        // Choose correct total cross-section array
        double sigma_tot = 0.0;
        if (projectile.name == "electron")
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion")
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }

        
        double nu = sigma_tot * velocity;
        total_nu += nu;
    }

    return total_nu / num_particles;
}

void CollisionHandler::coll_rate(Species &projectile, Species &target_gas)
{
    vec<double> rate(domain.ni);

    for (const Particle &part : projectile.part_list)
    {
        // Projectile lab-frame velocity (normalized back to physical units)
        double v_projx = part.vx * domain.vel_norm;
        double v_projy = part.vy * domain.vel_norm;
        double v_projz = part.vz * domain.vel_norm;

        double velocity = 0.0;

        if (projectile.name == "electron" || projectile.name == "ebeam")
        {
            // Cold gas approximation: v_rel ≈ v_proj
            double v_sqr = v_projx*v_projx + v_projy*v_projy + v_projz*v_projz;
            velocity = sqrt(v_sqr);
        }
        else
        {
            // Sample neutral velocity (Maxwellian at target_gas.temp)
            double vxg = Init::SampleVel(target_gas, target_gas.temp);
            double vyg = Init::SampleVel(target_gas, target_gas.temp);
            double vzg = Init::SampleVel(target_gas, target_gas.temp);

            // Relative velocity
            double v_relx = v_projx - vxg;
            double v_rely = v_projy - vyg;
            double v_relz = v_projz - vzg;
            velocity = sqrt(v_relx*v_relx + v_rely*v_rely + v_relz*v_relz);
        }

        // Energy in eV using mass of projectile species
        double energy = (0.5 * projectile.defaultmass * velocity * velocity) / Const::eV;

        // Get energy index
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        // Choose cross section
        double sigma_tot = 0.0;
        if (projectile.name == "electron" || projectile.name == "ebeam")
        {
            sigma_tot = sigma_tot_e[energy_index];
        }
        else if (projectile.name == "ion")
        {
            sigma_tot = sigma_tot_pion[energy_index];
        }
        else if (projectile.name == "beam" || projectile.name == "negion") //negative ion beam and negative ion
        {
            sigma_tot = sigma_tot_nion[energy_index];
        }

        // Local collision rate (per-particle contribution)
        double coll_rate = sigma_tot * velocity * domain.DT;

        // Deposit onto spatial grid
        domain.Scatter(part.x, coll_rate, projectile.coll_rate);
    }
}
