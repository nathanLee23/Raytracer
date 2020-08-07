#pragma once
#include "Obj.h"

class Sphere : public Obj {
public:
	Vec3 c;
	float r;

	Sphere(Vec3 c, float r);
	float intersect(Ray ray);
	Vec3 normal(Vec3 p);
};

