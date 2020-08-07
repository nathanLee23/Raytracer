#include "Plane.h"
#include "constants.h"

Plane::Plane(Vec3 a, Vec3 _n) {
	n = _n.normalized();
	k = n.dot(a);
}

float Plane::intersect(Ray ray) {
	return (k - n.dot(ray.o)) / n.dot(ray.d);
}

Vec3 Plane::normal(Vec3 p) {
	return n;
}