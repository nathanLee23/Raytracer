#include "Plane.h"
#include "constants.h"

Plane::Plane(Vec3 a, Vec3 _n) {
	n = _n.normalized();
	k = n.dot(a);
}

float Plane::intersect(Ray ray) {
	float d = n.dot(ray.d);
	if (d == 0.0) {
		return -1.0;
	}
	return (k - n.dot(ray.o)) / d;
}

Vec3 Plane::normal(Vec3 p) {
	return n;
}