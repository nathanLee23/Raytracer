#pragma once
#include "Obj.h"

class Plane : public Obj {
public:
	double k; // k = a.dot(n)
	Vec3 n;

	Plane(Vec3 a, Vec3 n);
	double intersect(Ray ray);
	Vec3 normal(Vec3 p);
};

