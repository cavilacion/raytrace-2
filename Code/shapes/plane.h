#ifndef PLANE_H_
#define PLANE_H_

#include "../object.h"

class Plane: public Object
{
    public:
        Plane(Point const  &point, Vector normal);

		virtual Hit intersect(Ray const &ray);

		virtual Color mapPointToTextureCoordinates(Point p);

		Point const point;

		Vector const normal;
};

#endif
