#pragma once
#include "Vec3.h"

enum Surface {
	diffuse,
	specular,
	reflective
};

struct Material {
	Vec3 albedo;
	float emission;
	Surface surface;
};