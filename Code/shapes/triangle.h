#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../object.h"

class Triangle: public Object
{
    public:
        Triangle(Point const  &A, Point const &B, Point const &C);

        virtual Hit intersect(Ray const &ray);

        Point const A, B, C; // three vertices
        
        Vector const N;
};

#endif
