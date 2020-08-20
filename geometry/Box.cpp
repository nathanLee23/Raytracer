#include <algorithm>
#include <cmath>

#include "Box.h"
#include "constants.h"
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

#include<iostream>

Vec3 Box::normal(Vec3 p) {
	p = -p;
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