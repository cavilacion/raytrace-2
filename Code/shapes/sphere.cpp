#include "sphere.h"
#include <iostream>
#include <cmath>

using namespace std;

Hit Sphere::intersect(Ray const &ray)
{
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, position, r
    * Sought: intersects? if true: *t
    *
    * Insert calculation of ray/sphere intersection here.
    *
    * You have the sphere's center (C) and radius (r) as well as
    * the ray's origin (ray.O) and direction (ray.D).
    *
    * If the ray does not intersect the sphere, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/

    Vector OC = (ray.O-position);
    Vector v = ray.D;
		
	double t, a, b, c;
	a=1.0;
	b=2.0*v.dot(OC);
	c=OC.dot(OC)-r*r;
	double D = b*b - 4.0*a*c;
    if (D<0) {
        return Hit::NO_HIT();
    }

    /****************************************************
    * RT1.2: NORMAL CALCULATION
    *
    * Given: t, C, r
    * Sought: N
    *
    * Insert calculation of the sphere's normal at the intersection point.
    ****************************************************/
		
	double t0 = (-b - sqrt(D) ) / (2*a);
	double t1 = (-b + sqrt(D) ) / (2*a);
	t = (t0 < t1 && t0 >= 0) ? t0 : t1;
	if (t<=0) {
        return Hit::NO_HIT();
    }         
	Point P = ray.at(t);
    Vector N = (P - position).normalized();
    return Hit(t,N);
}

Sphere::Sphere(Point const &pos, double radius)
:
    position(pos),
    r(radius)
{}
