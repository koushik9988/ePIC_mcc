#ifndef _SOLVERS_
#define _SOLVERS_

#include <vector>
#include <stdexcept>
#include "linalg.h"
#include "function.h"

template<typename T>
class Matrix;

template<typename T>
class vec;



// Function declarations
template<typename T>
vec<T> cg(Matrix<T> &A, vec<T> &x, vec<T> &b, int max_iteration, double tolerance);

template<typename T>
vec<T> pcg(Matrix<T> &A, vec<T> &x, vec<T> &b, int max_iteration, double tolerance);

template<typename T>
vec<T> direct(vec<T> rho, int n);

template<typename T>
vec<T> gausselimination(Matrix<T> &A, vec<T> &b);

//==========

// Conjugate Gradient method
template<typename T>
vec<T> cg(Matrix<T> &A, vec<T> &x, vec<T> &b, int max_iteration, double tolerance)
{
    if(symmetrycheck(A))
    {
        vec<T> r0 = b - A * x;
        vec<T> p0 = r0;

        for (int i = 0; i < std::min(b.getSize(), max_iteration); ++i)
        {
            if (r0.norm() < tolerance)
            {
                break;
            }

            vec<T> Ap0 = A * p0;
            double r02 = r0 * r0;
            double pAp = p0 * Ap0;
            double alpha = r02 / pAp;
            vec<T> x1 = x + p0 * alpha;
            vec<T> r1 = r0 - Ap0 * alpha;
            double r12 = r1 * r1;
            double beta = r12 / r02;
            vec<T> p1 = r1 + p0 * beta;

            x = x1;
            r0 = r1;
            p0 = p1;
        }
    }
    else
    {
        throw std::invalid_argument("fatal error! the Matrix is not symmetric,pcg/cg solver requre a symmetric matrix ");
    }

    return x;
}

// Preconditioned Conjugate Gradient method
template<typename T>
vec<T> pcg(Matrix<T> &A, vec<T> &x, vec<T> &b, int max_iteration, double tolerance)
{
    vec<T> r = b - A * x;
    Matrix<T> M_inv = diagonalPreconditioner(A);
    vec<T> w = M_inv * r;
    vec<T> v = w;
    double alpha = w * r;

    for (int iter = 0; iter < std::min(A.row * A.col, max_iteration); ++iter)
    {
        if (r.norm() < tolerance)
        {
            break;
        }

        vec<T> u = A * v;
        double t = alpha / (v * u);
        x = x + v * t;
        r = r - u * t;
        w = M_inv * r;
        double beta = w * r;
        double s = beta / alpha;
        v = w + v * s;
        alpha = beta;
    }

    return x;
}

// Direct solution using tridiagonal matrix solver
template<typename T>
vec<T> direct(vec<T> rho, int n)
{
    vec<T> a(n), b(n), c(n), x(n);

    for (int i = 1; i < n - 1; ++i)
    {
        a(i) = 1;
        b(i) = -2;
        c(i) = 1;
    }

    a(0) = 0;
    b(0) = 1;
    c(0) = 0;
    a(n - 1) = 0;
    b(n - 1) = 1;
    c(n - 1) = 0;

    x(0) = rho(0);
    x(n - 1) = rho(n - 1);

    for (int i = 1; i < n - 1; ++i)
    {
        x(i) = rho(i);
    }

    c(0) /= b(0);
    x(0) /= b(0);

    for (int i = 1; i < n; ++i)
    {
        double id = (b(i) - c(i - 1) * a(i));
        c(i) /= id;
        x(i) = (x(i) - x(i - 1) * a(i)) / id;
    }

    for (int i = n - 2; i >= 0; --i)
    {
        x(i) = x(i) - c(i) * x(i + 1);
    }

    return x;
}


/*
template<typename T>
vec<T> gausselimination(Matrix<T> &A, vec<T> &x, vec<T> &b)
//Matrix<T> gausselimination(Matrix<T> &A, vec<T> &x, vec<T> &b)
{
    int row = A.row;
    int col = A.col;
    double c;
    Matrix<T> M(row,col);
    for(int i = row-1 ; i > 1; i--)
    {
        for(int j = 0 ; j < col; j++)
        {
            c = (A(i,j)/A(i-1,j));
            if(i>j)
            {
                A(i,j) = A(i,j) - A(i-1,j)*c;
                b(i) = b(i) - b(i-1)*c;
            }
        }
    }

    for(int i = 0 ; i < row; i++)
    {
        for(int j = 0 ; j < col; j++)
        {
            {
                A(i,j) = A(i,j)/A(i,i);
                b(i) = b(i)/A(i,i);
            }
        }
    }

    // Back Substitution
    for(int i = row-1; i >= 0; i--)
    {
        x(i) = b(i);
        for(int j = i+1; j < row; j++)
        {
            x(i) -= A(i,j) * x(j);
        }
    }
    
    return x;

}
*/



template<typename T>
vec<T> gausselimination(Matrix<T> &A,vec<T> &b)
{
    int row = A.row;
    int col = A.col;
    double c;

    // Forward Elimination
    for(int i = 0; i < row-1; i++)
    {
        for(int k = i+1; k < row; k++)
        {
            c = A(k,i) / A(i,i);
            for(int j = 0; j < col; j++)
            {
                A(k,j) -= c * A(i,j);
            }
            b(k) -= c * b(i);
        }
    }
    vec<T> x(b.getSize());
    // Back Substitution
    for(int i = row-1; i >= 0; i--)
    {
        x(i) = b(i);
        for(int j = i+1; j < row; j++)
        {
            x(i) -= A(i,j) * x(j);
        }
        //x(i) /= A(i,i);
    }

    return x;
}

#endif // _SOLVERS_
