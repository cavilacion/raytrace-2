#include "triangle.h"

#include <cmath>

#define eps 0.0001

Hit Triangle::intersect(Ray const &ray)
{
	// https://stackoverflow.com/questions/995445/determine-if-a-3d-point-is-within-a-triangle
		if (ray.D.dot(N) == 0) { // ray is parallel to triangle
			return Hit::NO_HIT();
		}
			
    double t = (A - ray.O).dot(N) / (ray.D.dot(N)) ; // distance!
    Point P = ray.O + t*ray.D; 											 // intersection!
    Vector N1 = (B-A).cross(P-A).normalized();
    Vector N2 = (C-B).cross(P-B).normalized();
    Vector N3 = (A-C).cross(P-C).normalized();
    if (fabs(1.0-N1.dot(N2)) > eps ) {
			return Hit::NO_HIT();
		} 
		if (fabs(1.0-N2.dot(N3)) > eps ) {
			return Hit::NO_HIT();
		} 

    return Hit(t, N);
}

Triangle::Triangle(Point const &A, Point const &B, Point const &C)
:
		A(A),
		B(B),
		C(C),
		N((B - A).cross(C - A).normalized())
{}
