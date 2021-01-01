#pragma once

#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Vec3.h"
#include "Ray.h"
#include "Material.h"
#include "globals.h"


float sgn(float x) {
	return (float)((x > 0) - (x < 0));
}


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
//
//std::uniform_real_distribution<float> eta(0, 1);
//std::uniform_int_distribution<int> side(0, 5);
class Box : public Obj {
public:
	Vec3 min;
	Vec3 max;

	Box(Vec3 min, Vec3 max) : min(min), max(max) {
	}

	float intersect(Ray& ray) {
		float t1 = (min.x - ray.o.x) / ray.d.x;
		float t2 = (max.x - ray.o.x) / ray.d.x;

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

	Vec3 normal(Vec3 p) {
		Vec3 c = (max + min) / 2.0f;
		Vec3 cToP = p - c;

		// Can optimize this further
		int argMax = 0;
		for (int i = 1; i < 3; i++) {
			if (abs(cToP[i])*(max[argMax] - min[argMax]) > abs(cToP[argMax]) * (max[i] - min[i])) {
				argMax = i;
			}
		}

		return Vec3(
			(argMax == 0) * sgn(cToP.x),
			(argMax == 1) * sgn(cToP.y),
			(argMax == 2) * sgn(cToP.z)
		);
	}


	Vec3 samplePoint(std::mt19937 gen) {

		std::uniform_real_distribution<float> rx(min.x, max.x);
		std::uniform_real_distribution<float> rz(min.z, max.z);

		return Vec3(rx(gen), max.y, rz(gen));
	}
};

class Triangle : public Obj {
public:
	Vec3 a;
	Vec3 b;
	Vec3 c;
	Vec3 n;

	Triangle(Vec3 a, Vec3 b, Vec3 c) : a(a), b(b), c(c), n((b - a).cross(c - b).normalized()) {
	}

	float intersect(Ray& ray) {
		return 0.0;
	}
	Vec3 normal(Vec3 p) {
		return n;
	}
	Vec3 samplePoint(std::mt19937 gen) {
		return Vec3();
	}
};

class Mesh : public Obj {
public:
	std::vector<Triangle> triangles;
	Box bounds = Box(Vec3(), Vec3());

	Mesh() {
	}

	void  add_triangle(Triangle t) {
		triangles.push_back(t);
	}
	float intersect(Ray& ray) {
		if (triangles.empty()) {
			return INFINITY;
		}
		float min_t = INFINITY;
		Triangle min_triangle = triangles.front();
		for (Triangle triangle : triangles) {
			float t = triangle.intersect(ray);
			if (t < min_t) {
				min_t = t;
				min_triangle = triangle;
			}
		}
		return min_t;
	}
	Vec3 normal(Vec3 p) {
		return Vec3();
	}
	Vec3 samplePoint(std::mt19937 gen) {
		return Vec3();
	}
};

int filter_extra(std::istringstream &s) {
	std::string x;
	if (s >> x) {
		return std::stoi(x.substr(0, x.find_first_of('/')));
	} else {
		return -1;
	}

}

Mesh load_mesh(std::string file_name) {
	std::ifstream ifile(file_name);
	if (errno) {
		throw std::system_error(errno, std::system_category(), "Failed to open: " + file_name);
	} else {
		std::cout << "Successfully opened file " + file_name << std::endl;
	}
	Mesh mesh = Mesh();

	
	std::vector<Vec3> vertices;
	while (!ifile.eof()) {
		if (isspace(ifile.peek())) {
			ifile.get();
			continue;
		}
		std::string line;
		std::getline(ifile, line);
		if (line.size() < 2) { // Too short
		} else if (line.rfind("#", 0) == 0) { // comment
		} else if (line.rfind("v ", 0) == 0) { // vertex
			Vec3 vertex;
			sscanf_s(line.c_str(), "v %f %f %f", &vertex.x, &vertex.y, &vertex.z);
			vertices.push_back(vertex);
		} else if (line.rfind("f ", 0) == 0) {//face
			std::istringstream line_ss(line);
			std::string word;
			line_ss >> word;

			int a = filter_extra(line_ss);
			int b = filter_extra(line_ss);
			int c = filter_extra(line_ss);

			mesh.add_triangle(Triangle(vertices.at(a - 1), vertices.at(b - 1), vertices.at(c - 1)));
			b = c;
			while (c = filter_extra(line_ss) >= 0) {
				mesh.add_triangle(Triangle(vertices.at(a - 1), vertices.at(b - 1), vertices.at(c - 1)));
				b = c;
			}
		}
	}
	return mesh;
}

class Plane : public Obj {
public:
	float k; // k = a.dot(n)
	Vec3 n;

	Plane(Vec3 a, Vec3 n_) {
		n = n_.normalized();
		k = n.dot(a);
	}

	float intersect(Ray& ray) {
		return (k - n.dot(ray.o)) / n.dot(ray.d);
	}

	Vec3 normal(Vec3 p) {
		return n;
	
	}
	Vec3 samplePoint(std::mt19937 gen) {
		return Vec3();
	}
};

class Sphere : public Obj {
public:
	Vec3 c;
	float r;

	Sphere(Vec3 c, float r) : c(c), r(r) {
	}
	float intersect(Ray& ray) {
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
	Vec3 normal(Vec3 p) {
		return (p - c) / r;
	}
	Vec3 samplePoint(std::mt19937 gen) {
		return Vec3();
	}
};