#include "Ray.h"
#include "constants.h"
#include "Vec3.h"
#include "Sphere.h"


Sphere::Sphere(Vec3 c, float r) : c(c), r(r) {
}

float Sphere::intersect(Ray ray) {
	float B = 2 * ray.d.dot(ray.o-c);
	float C = (ray.o-c).sqrNorm() - r * r;

	float sqrDiscr = B*B - 4 * C;
	if (sqrDiscr < 0.0f) {
		return INFINITY;
	}
	float discr = sqrt(sqrDiscr);
	float t1 = -B - discr;
	float t2 = -B + discr;
	return t1 > 2.0f*EPS ? t1 / 2.0f : (t2 > 2.0f*EPS ? t2 / 2.0f : INFINITY);
}

Vec3 Sphere::normal(Vec3 p) {
	return (p - c)/r;
}
