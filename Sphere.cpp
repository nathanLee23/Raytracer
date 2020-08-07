#include "Sphere.h"

#include "Ray.h"
#include "constants.h"
#include <algorithm>
#include "Vec3.h"

Sphere::Sphere(Vec3 c, float r) : c(c), r(r) {
}

float Sphere::intersect(Ray ray) {
	float B = 2 * ray.d.dot(ray.o-c);
	float C = (ray.o-c).sqrNorm() - r * r;

	float sqrDiscr = B*B - 4 * C;
	if (sqrDiscr < 0.0) {
		return -1.0;
	}
	float discr = sqrt(sqrDiscr);
	float t1 = -B - discr;
	float t2 = -B + discr;
	return t1 > EPS ? t1 / 2.0 : (t2 > EPS ? t2 / 2.0 : -1.0);	
}

Vec3 Sphere::normal(Vec3 p) {
	return (p - c)/r;
}
