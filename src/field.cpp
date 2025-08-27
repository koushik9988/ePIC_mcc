#include "field.h"

//mo problem here(numerical error problem )
//void FieldSolve::SolvePotDirect(double *x, double *rho)

void FieldSolve::Spectral()
{
    int ni = domain.ni;
    double L = domain.xL;
    int nr = (ni/2) + 1;
    double norm = 1.0 /ni;
 
    double *rho = fftw_alloc_real(ni);
    double *phi = fftw_alloc_real(ni);
   
    fftw_complex *rho_k= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nr));
    fftw_complex *phi_k= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nr));

    for(int i = 0 ; i < ni ; i++)
    {
        rho[i] =  -domain.rho(i)*((domain.L*domain.L)/(domain.LDe*domain.LDe));
    }

    // Forward FFT
    fftw_plan forward = fftw_plan_dft_r2c_1d(ni, rho , rho_k, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    //display::print(phi_k[0][0]);
    // Solve in Fourier space
    for (int k = 0; k <= ni/2; ++k)
    {
        //double kx = 2.0 * M_PI * k / L;
        //double denom = (k == 0) ? 1.0 : -(kx * kx); // Avoid divide-by-zero at k=0
        double dx = domain.dx;
        double kx = 2.0 * M_PI * k / L;
        double denom = (k == 0) ? 1.0 : -kx*kx;//- (4/(dx*dx))*sin(kx*dx*0.5)*sin(kx*dx*0.5);

        phi_k[k][0] = (k == 0) ? 0.0 : rho_k[k][0] / denom;
        phi_k[k][1] = (k == 0) ? 0.0 : rho_k[k][1] / denom;

        domain.dft_k(k) = kx;
    }

    //display::print(phi_k[0][0]);

    for( int i = 0 ; i < nr ; i++)
    {
        domain.dft_value(i) = sqrt(rho_k[i][0]*rho_k[i][0] + rho_k[i][1]*rho_k[i][1]);
    }

    // Inverse FFT
    fftw_plan backward = fftw_plan_dft_c2r_1d(ni, phi_k, phi, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    domain.phi = 0;
    for (int i = 0; i < ni; ++i)
    {
        domain.phi(i) = phi[i] * norm;
    }

    // Cleanup
    fftw_free(rho_k);
    fftw_free(phi_k);
    fftw_free(rho);
    fftw_cleanup();
    
}

void FieldSolve::Direct(int ts)
{
    /* Set coefficients, precompute them*/
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> rhov(ni);

    double frequency = 13.56e6;
    double omega = 2.0 * Const::PI * frequency;

    double volatge = 100.0;  // Amplitude of the electric field pulse
    
    if(domain.bc == "pbc")
    {  
        rhov(0) = 0;//domain.vL;// 10*sin(omega*ts*domain.DT*1)*(Const::eV/(Const::K_b*Const::EV_to_K));
        rhov(ni-1) = 0;//100*cos(omega*ts*domain.DT*10);   
    }
    if(domain.bc == "open")
    {
        //rhov(0) = domain.vL * domain.phi_norm;
        rhov(0) = domain.vL*domain.phi_norm ;//* sin(omega*ts*domain.DT); 
        rhov(ni-1) = domain.vR * domain.phi_norm;
    }
           
    for( int i = 1 ; i < ni-1 ; i++)
    {
        rhov(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
        //rhov(i) += 0.1*sin((2*3.14*i)/domain.xL)*((domain.L*domain.L)/(domain.LD*domain.LD));
    }

    vec<double> sol = direct(rhov,ni);

    domain.phi = 0;
    domain.phi = sol;
}


void FieldSolve::GaussElim()
{
    double L2;
    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    Matrix<double> A(n,n);

   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
    }
    
    vec<double> sol = gausselimination(A,b);

    domain.phi = 0;
    
    domain.phi = sol;
}


void FieldSolve::pcgsolver()
{
    double L2;
    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    vec<double> x(n);
    Matrix<double> A(n,n);

    std::vector<bool> fixed_node(n, false); 
    // Set boundary conditions
   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
        x(i) = domain.phi(i);
    }

    Matrix<double> A0 = slice(A,1,n-1,1,n-1);
    vec<double> b0 = slice(b,1,n-1);
    vec<double> x0 = slice(x,1,n-1);

    //symmetrycheck(A0);

    vec<double> x_sol = cg(A0,x0,b0,domain.max_iteration, domain.tolerance);
   
    vec<double> v_norm = A0*x_sol - b0;

    //reset phi
    domain.phi = 0;

    domain.phi(0) = b(0);
    domain.phi(n-1) = b(n-1);
    
    domain.phi = x_sol;

    for(int i = 1; i< n-1; i++)
    {
        domain.phi(i) = x_sol(i-1);
    }
}


void FieldSolve::cgsolver()
{
    double L2;

    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    vec<double> x(n);
    Matrix<double> A(n,n);

    std::vector<bool> fixed_node(n, false); 
    // Set boundary conditions
   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));;
        x(i) = domain.phi(i);
    }

    Matrix<double> A0 = slice(A,1,n-1,1,n-1);
    vec<double> b0 = slice(b,1,n-1);
    vec<double> x0 = slice(x,1,n-1);

    //symmetrycheck(A0);

    vec<double> x_sol = cg(A0,x0,b0,domain.max_iteration, domain.tolerance);
   
    vec<double> v_norm = A0*x_sol - b0;

    //display::print(x_sol.getSize());

    domain.norm = v_norm.norm();

    domain.phi(0) = b(0);
    domain.phi(n-1) = b(n-1);
    
    for(int i = 1; i< n-1; i++)
    {
        domain.phi(i) = x_sol(i-1);
    }
}

//void FieldSolve::CalculateEfield(double *phi, double *ef)
void FieldSolve::CalculateEfield()
{
    
        domain.ef = 0;
        for(int i=1; i<domain.ni-1; i++)
        {
            domain.ef(i) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(i+1)-domain.phi(i-1))/(2*domain.dx);
        }
        
        /*for continous bounndary
        the point 0 and ni-1 is same */
        if(domain.bc =="pbc")
        {
            domain.ef(0) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(1)-domain.phi(domain.ni-2))/(2*domain.dx);
            domain.ef(domain.ni-1) = domain.ef(0);
        }
        else if(domain.bc == "open" || domain.bc == "rbc")
        {
            domain.ef(0) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(1)- domain.phi(0))/(domain.dx);
            domain.ef(domain.ni-1) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(domain.ni-1)-domain.phi(domain.ni-2))/(domain.dx);
        }
}


void FieldSolve::AddPerturbation(int ts, int mode, double perturb_amplitude)
{
    double k0 = 2.0 * M_PI * mode / domain.xL;  
    double amplitude = perturb_amplitude;
    double phase = 0.0;

    for (int i = 0; i < domain.ni; i++)
    {
        double x = i * domain.dx;
        domain.phi(i) += amplitude * cos(k0 * x + phase);
    }
}



void FieldSolve::PotentialSolver(int ts)
{
    if(domain.SolverType == "direct") Direct(ts);
    if(domain.SolverType == "cg") cgsolver();
    if(domain.SolverType == "pcg") pcgsolver();
    if(domain.SolverType == "spectral") Spectral();
}
