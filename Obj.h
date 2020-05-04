#pragma once
#include "Vec3.h"
#include "Ray.h"

class Obj {
public:
	Vec3 albedo = Vec3(255.0);
	double emission = 0.0;

	// Returns 
	virtual double intersect(Ray ray) = 0; 
	virtual Vec3 normal(Vec3 p) = 0;
};

