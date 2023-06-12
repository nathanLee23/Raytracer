#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include "globals.h"
#include "sampler.h"

// Also called pinhole camera
class PerspectiveCamera {
public:
	float fov;
	PerspectiveCamera(float fov) : fov(fov) {
	}

	// Transform pixel coordinates to perspective rays
	// Camera is at (0,0,0) facing (0,0,-1)
	void pixelToCameraRay(float px, float py, Ray& ray, Pcg& gen) const {
		float x = (2.0f * px - WIDTH) / WIDTH * tanf(fov * M_PI / 180.0f / 2.0f);
		float y = (2.0f * py - HEIGHT) / HEIGHT * tanf(((float)HEIGHT) / WIDTH * fov * M_PI / 180.0f / 2.0f);
		float z = -1.0f;

		ray.o = Vec3(0.0f, 1.0f, 3.0f);
		ray.d = Vec3(x, -y, z).normalized();
	}
};

Vec3 uniformSampleUnitDisk(Pcg& gen) {
	float r = gen.Uniform();
	float theta = gen.Uniform() * 2 * M_PI;

	return Vec3(r*std::cos(theta), r*std::sin(theta), 0.0f);
}

class ThinLensCamera {
public:
	float focal_length = 2.4f;
	float lens_radius = 0.09f;
	float fov = 60.0f;

	ThinLensCamera() {
	}

	void pixelToCameraRay(float px, float py, Ray& ray, Pcg& gen) const {
		float x = (2.0f * px - WIDTH) / WIDTH * tanf(fov * M_PI / 180.0f / 2.0f);
		float y = (2.0f * py - HEIGHT) / HEIGHT * tanf(((float)HEIGHT) / WIDTH * fov * M_PI / 180.0f / 2.0f);
		float z = -1.0f;

		ray.o = Vec3(0.0f, 1.0f, 3.0f);
		ray.d = Vec3(x, -y, z);

		float t = -focal_length / ray.d.z;
		Vec3 focal_point = ray.o + t * ray.d;

		ray.o += uniformSampleUnitDisk(gen) * lens_radius;
		ray.d = (focal_point - ray.o).normalized();
	}
};


class OrthogonalCamera {
public:

	float stretch_x;
	float stretch_y;
	OrthogonalCamera(float stretch_x, float stretch_y) : stretch_x(stretch_x), stretch_y(stretch_y) {
	}
	OrthogonalCamera(float stretch) : OrthogonalCamera(stretch, stretch) {
	}

	// Transform pixel coordinates to perspective rays
	// Camera is at (0,0,0) facing (0,0,-1)
	void pixelToCameraRay(float px, float py, Ray& ray, Pcg& gen) const {
		// Convert to NDC and stretch
		ray.o = Vec3(0.0f, 1.0f, 3.0f) +
			Vec3((2.0f * px - WIDTH) / WIDTH * stretch_x, -(2.0f * py - HEIGHT) / HEIGHT * stretch_y, 0.f);

		ray.d = Vec3(0.0f, 0.0f, -1.0f);
	}
};