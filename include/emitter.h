#include "species.h"
#include "domain.h"
#include "init.h"

class Domain;
class Species;

struct EmitterParams
{
    double emitter_loc;
    double temp;
    int numparticle;
    double vdrift;
    int species_id1;
    int species_id2;
};

class Emitter
{
    public:
    Emitter(const EmitterParams &params, Domain &domain);
    void inject(std::vector<Species> &species_list); 

    private:
    EmitterParams params;
    //std::vector<Species> species_list;
    Domain &domain;
};
