#include "plane.h"

#include <cmath>

Hit Plane::intersect(Ray const &ray)
{
	float denom = normal.dot(ray.D);
	if (fabs(denom) <= 0.0001f) // your favorite epsilon
	{
		return Hit::NO_HIT();
	}
	float t = (point - ray.O).dot(normal) / denom;
	
  return Hit(t, normal);
}

Plane::Plane(Point const  &point, Vector normal)
:
	point(point),
	normal(normal.normalized())
{}
