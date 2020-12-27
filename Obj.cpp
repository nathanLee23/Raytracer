#include <algorithm>
#include <random>
#include <cmath>

#include "Obj.h"
#include "globals.h"
#include "Vec3.h"

float sgn(float x) {
	return (x > 0) - (x < 0);
}

Box::Box(Vec3 min, Vec3 max) : min(min), max(max) {
}

float Box::intersect(Ray& ray) {
	float t1 = (min.x - ray.o.x)/ray.d.x;
	float t2 = (max.x - ray.o.x)/ray.d.x;
	
	float tmin = std::min(t1, t2);
	float tmax = std::max(t1, t2);

	t1 = (min.y - ray.o.y) / ray.d.y;
	t2 = (max.y - ray.o.y) / ray.d.y;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));
	t1 = (min.z - ray.o.z) / ray.d.z;
	t2 = (max.z - ray.o.z) / ray.d.z;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));
	//for (int i = 1; i < 3; i++) {
	//	t1 = (min[i] - ray.o[i])/ray.d[i];
	//	t2 = (max[i] - ray.o[i])/ray.d[i];

	//	tmin = std::max(tmin, std::min(t1, t2));
	//	tmax = std::min(tmax, std::max(t1, t2));
	//}
	return tmax >= tmin ? (tmin > EPS ? tmin : tmax) : INFINITY;
}

Vec3 Box::normal(Vec3 p) {
	Vec3 c = (max + min) / 2.0f;
	Vec3 cToP = p - c;

	// Can optimize this further
	int argMax = 0;
	for (int i = 1; i < 3; i++) {
		if (abs(cToP[i])*(max[argMax] - min[argMax]) > abs(cToP[argMax]) * (max[i] - min[i]) ) {
			argMax = i;
		}
	}

	return Vec3(
		(argMax == 0) * sgn(cToP.x),
		(argMax == 1) * sgn(cToP.y),
		(argMax == 2) * sgn(cToP.z)
	);
}
std::uniform_real_distribution<float> eta(0,1);
std::uniform_int_distribution<int> side(0, 5);
Vec3 Box::samplePoint(std::mt19937 gen) {

	std::uniform_real_distribution<float> rx(min.x, max.x);
	std::uniform_real_distribution<float> rz(min.z, max.z);
	
	return Vec3(rx(gen), max.y, rz(gen));
}
Vec3 Plane::samplePoint(std::mt19937 gen) {
	return Vec3();
}
Vec3 Sphere::samplePoint(std::mt19937 gen) {
	return Vec3();
}
Plane::Plane(Vec3 a, Vec3 _n) {
	n = _n.normalized();
	k = n.dot(a);
}

float Plane::intersect(Ray& ray) {
	return (k - n.dot(ray.o)) / n.dot(ray.d);
}

Vec3 Plane::normal(Vec3 p) {
	return n;
}

Sphere::Sphere(Vec3 c, float r) : c(c), r(r) {
}

float Sphere::intersect(Ray& ray) {
	float B = 2 * ray.d.dot(ray.o - c);
	float C = (ray.o - c).sqrNorm() - r * r;

	float sqrDiscr = B * B - 4 * C;
	if (sqrDiscr < 0.0f) {
		return INFINITY;
	}
	float discr = sqrt(sqrDiscr);
	float t1 = (-B - discr) / 2.0f;
	//float t2 = (-B + discr)/2.0f;
	return t1 > EPS ? t1 : (-B + discr) / 2.0f;
}

Vec3 Sphere::normal(Vec3 p) {
	return (p - c) / r;
}