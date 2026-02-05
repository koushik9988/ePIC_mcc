#ifndef _EXTRA_
#define _EXTRA_

#include <iostream>

inline bool string_to_bool(std::string pflag)
{
    bool val;
    if(pflag == "true" || pflag == "True" || pflag == "T" || pflag == "t" || pflag == "1" || pflag == "yes" || pflag == "Yes" || pflag == "Y" || pflag == "y")
    {
        val = true;
    }
    else
    {
        val = false;
    }
    return val;
}

#endif