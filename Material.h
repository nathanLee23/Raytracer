#pragma once
#include "Vec3.h"

enum Surface {
	diffuse,
	specular,
	reflective,
	varnish
};

struct Material {
	Vec3 albedo;
	float emission;
	Surface surface;
};


class _Material {
	virtual Vec3 brdf(Vec3 wi, Vec3 wo) = 0;

	virtual float sampleDirection(Vec3 wi) = 0;

	virtual float pdf(Vec3 wi, Vec3 wo) = 0;
};

class Diffuse : public _Material {
	Vec3 albedo;
	Vec3 emission;

	Vec3 brdf(Vec3 wi, Vec3 n, Vec3 wo) {
		return 1.0f / albedo;
	}

	float sampleDirection(Vec3 wi) {

	}

	float pdf(Vec3 wi, Vec3 wo) {
		//return ;
	}
};


union __Material {

};