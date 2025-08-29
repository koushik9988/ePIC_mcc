#include <iostream>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <fstream>
#include "slap.h"
#ifdef ENABLE_PLOTTING
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

void SpectralSolver(vec<double> &rho_in, vec<double> &phi_out, int n, double L)
{
    int nr = (n/2) + 1;
    double norm = 1.0 / n;
    double dx = L/(n-1);
 
    double *rho = fftw_alloc_real(n);
    double *phi = fftw_alloc_real(n);
    fftw_complex *rho_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nr));
    fftw_complex *phi_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nr));

    for(int i = 0 ; i < n ; i++)
    {
        rho[i] = -rho_in(i);
    }

    // Forward FFT
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, rho , rho_k, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // Solve in Fourier space
    for (int k = 0; k <= n/2; ++k)
    {
        double kx = 2.0 * M_PI * k / L;
        double denom = (k == 0) ? 1.0 : - (4/(dx*dx))*sin(kx*dx*0.5)*sin(kx*dx*0.5); // Avoid divide-by-zero at k=0

        phi_k[k][0] = (k == 0) ? 0.0 : rho_k[k][0] / denom;
        phi_k[k][1] = (k == 0) ? 0.0 : rho_k[k][1] / denom;
    }

    // Inverse FFT
    fftw_plan backward = fftw_plan_dft_c2r_1d(n, phi_k, phi, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize
    for (int i = 0; i < n; ++i)
    {
        phi[i] *= norm;
    }

    for(int i = 0 ; i < n ; i++)
    {
        phi_out(i) = phi[i];
    }

    // Cleanup
    fftw_free(rho_k);
    fftw_free(phi_k);
    fftw_free(rho);
    fftw_cleanup();
}

int main()
{
    int n = 1024;
    double L = 100.0;
    double dx = L / (n-1);

    vec<double> x(n), rho(n), phi(n), phi_exact(n);

    for (int i = 0; i < n; ++i)
    {
        x(i) = i * dx;
        rho(i) = std::sin(2.0 * M_PI * x(i) / L) + std::cos(4.0 * M_PI * x(i) / L);
        phi_exact(i) = (L*L)/(4*M_PI*M_PI) * std::sin(2.0 * M_PI * x(i)/L) + (L*L)/(16*M_PI*M_PI) * std::cos(4.0 * M_PI * x(i) / L);
    }

    SpectralSolver(rho, phi, n, L);

#ifdef ENABLE_PLOTTING
    std::vector<double> x1(n);
    std::vector<double> phi_spectral(n);
    std::vector<double> phi_exact_sol(n);

    for (int i = 0; i < n; ++i)
    {
        x1[i] = x(i);
        phi_spectral[i] = phi(i);
        phi_exact_sol[i] = phi_exact(i);
    }
    
    plt::named_plot("spectral", x1, phi_spectral, "r-");
    plt::named_plot("exact", x1, phi_exact_sol, "b-");
    plt::legend();
    plt::xlabel("x");
    plt::ylabel("$phi$");
    plt::show();
#endif

#ifndef ENABLE_PLOTTING
    // Optional: Output results to console or file when plotting is disabled
    std::ofstream out("spectral_solver_results.txt");
    out << "x\tphi_spectral\tphi_exact\n";
    for (int i = 0; i < n; ++i)
    {
        out << x(i) << "\t" << phi(i) << "\t" << phi_exact(i) << "\n";
    }
    out.close();
    std::cout << "Spectral solver results written to spectral_solver_results.txt" << std::endl;
#endif

    return 0;
}