#pragma once
#include "Obj.h"

class Sphere : public Obj {
public:
	Vec3 c;
	double r;

	Sphere(Vec3 c, double r);
	double intersect(Ray ray);
	Vec3 normal(Vec3 p);
};

