#pragma once
#include "Obj.h"

class Plane : public Obj {
public:
	float k; // k = a.dot(n)
	Vec3 n;

	Plane(Vec3 a, Vec3 n);
	float intersect(Ray ray);
	Vec3 normal(Vec3 p);
};

