#include "emitter.h"

Emitter::Emitter(const EmitterParams &params, Domain &domain) : params(params), domain(domain) {}

void Emitter::inject(std::vector<Species> &species_list)
{
    for (int p = 0; p < params.numparticle; p++)
    {
        // In 1D, position is fixed at emitter_loc
        double x = params.emitter_loc;

        // First particle
        double vx1 = Init::SampleVel(species_list[params.species_id1], params.temp) + params.vdrift * domain.vel_norm;

        vx1 /= domain.vel_norm;

        species_list[params.species_id1].AddParticle(Particle(x, vx1, 0, 0,0));

        // If the two species indices are different, inject the second particle
        if (params.species_id1 != params.species_id2)
        {
            double vx2 = Init::SampleVel(species_list[params.species_id2], params.temp) + params.vdrift * domain.vel_norm;
            
            vx2 /= domain.vel_norm;

            species_list[params.species_id2].AddParticle(Particle(x, vx2, 0, 0,0));
        }
    }
}