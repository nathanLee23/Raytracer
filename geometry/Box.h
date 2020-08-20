#pragma once
#include "Obj.h"

class Box : public Obj {
public:
	Vec3 min;
	Vec3 max;

	Box(Vec3 min, Vec3 max);
	float intersect(Ray& ray);
	Vec3 normal(Vec3 p);
};

