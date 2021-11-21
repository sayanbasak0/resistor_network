#include "covernout.h"
#include <iostream>
#include <cmath>
#include <float.h>

int main() {
    std::cout << DBL_MAX << "\n"
              << 1/DBL_MAX << "\n"
              << DBL_MAX+(DBL_MAX/1000.0) << "\n"
              << ((DBL_MAX+(DBL_MAX/1000.0))==0) << "\n"
              << (DBL_MAX+DBL_MAX)*0 << "\n"
              << DBL_MAX/(DBL_MAX/10.0) << "\n"
              << (0.1/(1+0.1))*DBL_MAX << "\n"
              << (DBL_MAX/2)+(DBL_MAX/2)+(DBL_MAX/2) << "\n"
              << (DBL_MAX/2)+(DBL_MAX/1.5) << "\n"
              << sqrt(DBL_MAX)*sqrt(DBL_MAX)/DBL_MAX << "\n"
              << 1.2*DBL_MAX << "\n"
              << 10*DBL_MAX << "\n"
              << fmin(1.0,2.0) << "\n"
              << fmax(1,2) << "\n"
              << std::isfinite(1.2) << "\n"
              << std::isinf(1.2) << "\n"
              ;
}

// a*b/(a+b)
// a>b
// r = b/a
// r<1 
// r*a*a/((r+1)*a)
// r/(r+1)*a