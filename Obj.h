#pragma once
#include <random>

#include "Vec3.h"
#include "Ray.h"
#include "Material.h"

class Obj {
public:
	Material material;

	// Returns the t for which r.o+t*r.d intersects the object
	// Returns a negative number if no intersection exists
	virtual float intersect(Ray& ray) = 0;

	// If p is a point on the object normal() returns the normalized normal of the object at the point
	virtual Vec3 normal(Vec3 p) = 0;
	virtual Vec3 samplePoint(std::mt19937 gen) = 0;
};

class Box : public Obj {
public:
	Vec3 min;
	Vec3 max;

	Box(Vec3 min, Vec3 max);
	float intersect(Ray& ray);
	Vec3 normal(Vec3 p);
	Vec3 samplePoint(std::mt19937 gen);
};

class Plane : public Obj {
public:
	float k; // k = a.dot(n)
	Vec3 n;

	Plane(Vec3 a, Vec3 n);
	float intersect(Ray& ray);
	Vec3 normal(Vec3 p);
	Vec3 samplePoint(std::mt19937 gen);
};

class Sphere : public Obj {
public:
	Vec3 c;
	float r;

	Sphere(Vec3 c, float r);
	float intersect(Ray& ray);
	Vec3 normal(Vec3 p);
	Vec3 samplePoint(std::mt19937 gen);
};