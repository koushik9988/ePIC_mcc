#ifndef _MATH_
#define _MATH_

#include "slap.h"

//trapezoidal rule
inline double trape(vec<double> &field, double h, int n)
{
    double sum = 0.5 * (field(0) + field(n-1)); 

 
    for (int i = 1; i < n-1; ++i)
    {
        sum += field(i);
    }

    return h * sum;
}

//Simpson's 3/8 rule
inline double simpson_38(vec<double> &field, double h, int n)
{
    double sum = field(0) + field(n-1);

    for (int i = 1; i < n-1; ++i)
    {
        if (i % 3 == 0)
        {  
            sum += 2 * field(i);
        }
        else
        {
            sum += 3 * field(i);
        }
    }

    return (3 * h / 8) * sum;
}

//Simpson's 1/3 rule
inline double simpson_13(vec<double> &field, double h, int n)
{
    double sum = field(0) + field(n-1);

    for (int i = 1; i < n-1; ++i)
    {
        if (i % 2 == 0)
        { 
            sum += 2 * field(i);
        }
        else
        {
            sum += 4 * field(i);
        }
    }

    return (h / 3) * sum;
}

#endif
