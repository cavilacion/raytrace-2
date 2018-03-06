#ifndef TORUS_H_
#define TORUS_H_

#include "../object.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <float.h>

#define sqr(a) (a)*(a)
#define cub(a) (a)*sqr(a)
#define qrt(a) (a)*cub(a)
#define infty DBL_MAX

class Torus: public Object
{
    public:
        Torus(Point const &C, double r, double R);

        virtual Hit intersect(Ray const &ray);
        
        Point const C;
        
        double const R, r;
};

std::complex<double> cuberoot(std::complex<double> z);
std::complex<double> *solveQuadratic (const std::complex<double> *a);
std::complex<double> *solveQuarticComplex (const std::complex<double> *a);
int numRealSolutions (const double *a);
int solveQuarticReal (double &t, const double *c);
#endif
