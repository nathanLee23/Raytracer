#pragma once

#include "Vec3.h"

class Ray {
public:
	Vec3 o, d;
	Ray(Vec3 o, Vec3 d) : o(o), d(d.normalized()) {
	};
};
