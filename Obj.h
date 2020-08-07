#pragma once
#include "Vec3.h"
#include "Ray.h"
#include "Material.h"

class Obj {
public:
	Material material;

	// Returns the t for which r.o+t*r.d intersects the object
	// Returns a negative number if no intersection exists
	virtual float intersect(Ray ray) = 0;

	// If p is a point on the object normal() returns the normalized normal of the object at the point
	virtual Vec3 normal(Vec3 p) = 0;
};

