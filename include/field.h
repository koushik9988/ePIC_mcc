#ifndef _FIELD_H_
#define _FIELD_H_

#include <cstring>
#include "domain.h"
#include "linalg.h"
#include "function.h"
#include "solvers.h"
#include <fftw3.h>

class Domain; // Forward declaration of Domain class

class FieldSolve
{
public:
    // Constructor
    FieldSolve(Domain &domain):domain(domain){};

    enum class SolverType
    {
        DIRECT,
        PCG,
        CG,
        GS
    };

    void PotentialSolver(int ts);
    void CalculateEfield();
    void Direct(int ts);
    void pcgsolver();
    void cgsolver();
    void GaussElim();
    void Spectral();

private:
    Domain &domain;
};

#endif // _FIELD_H_
