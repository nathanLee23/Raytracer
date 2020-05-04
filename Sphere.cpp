#include "Sphere.h"

#include "Ray.h"
#include <algorithm>
#include "Vec3.h"

Sphere::Sphere(Vec3 c, double r) : c(c), r(r) {
}

double Sphere::intersect(Ray ray) {
	double B = 2 * ray.d.dot(ray.o-c);
	double C = (ray.o-c).sqrNorm() - r * r;

	double sqrDiscr = B*B - 4 * C;
	if (sqrDiscr > 0) {
	} else {
		return -1.0;
	}
	double discr = sqrt(sqrDiscr);
	double t1 = (-B - discr);
	double t2 = (-B + discr);
	// TODO: replace 0 with EPS
	return t1 > 0 ? t1 / 2.0 : (t2 > 0 ? t2 / 2.0 : -1.0);	
}

Vec3 Sphere::normal(Vec3 p) {
	return (p - c)/r;
}
