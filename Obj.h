#pragma once

#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

#include <embree3/rtcore.h>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include "Vec3.h"
#include "Ray.h"
#include "Material.h"
#include "globals.h"

#include "Sampler.h"


float sgn(float x) {
	return (float)((x > 0) - (x < 0));
}

class Obj {
public:
	// Returns the t for which r.o+t*r.d intersects the object
	// Returns a negative number if no intersection exists
	virtual float intersect(Ray& ray) = 0 ;

	virtual float samplePoint(Pcg& gen, Vec3& point_ret) const = 0;
};
//
// 
//class Box : public Obj {
//public:
//	Vec3 min;
//	Vec3 max;
//
//	Box(Vec3 min, Vec3 max) : min(min), max(max) {
//	}
//
//	float intersect(Ray& ray) {
//		float t1 = (min.x - ray.o.x) / ray.d.x;
//		float t2 = (max.x - ray.o.x) / ray.d.x;
//
//		float tmin = std::min(t1, t2);
//		float tmax = std::max(t1, t2);
//
//		t1 = (min.y - ray.o.y) / ray.d.y;
//		t2 = (max.y - ray.o.y) / ray.d.y;
//
//		tmin = std::max(tmin, std::min(t1, t2));
//		tmax = std::min(tmax, std::max(t1, t2));
//		t1 = (min.z - ray.o.z) / ray.d.z;
//		t2 = (max.z - ray.o.z) / ray.d.z;
//
//		tmin = std::max(tmin, std::min(t1, t2));
//		tmax = std::min(tmax, std::max(t1, t2));
//		//for (int i = 1; i < 3; i++) {
//		//	t1 = (min[i] - ray.o[i])/ray.d[i];
//		//	t2 = (max[i] - ray.o[i])/ray.d[i];
//
//		//	tmin = std::max(tmin, std::min(t1, t2));
//		//	tmax = std::min(tmax, std::max(t1, t2));
//		//}
//		return tmax >= tmin ? (tmin > EPS ? tmin : tmax) : INFINITY;
//	}
//
//	Vec3 normal(Vec3 p) {
//		Vec3 c = (max + min) / 2.0f;
//		Vec3 cToP = p - c;
//
//		// Can optimize this further
//		int argMax = 0;
//		for (int i = 1; i < 3; i++) {
//			if (abs(cToP[i])*(max[argMax] - min[argMax]) > abs(cToP[argMax]) * (max[i] - min[i])) {
//				argMax = i;
//			}
//		}
//
//		return Vec3(
//			(argMax == 0) * sgn(cToP.x),
//			(argMax == 1) * sgn(cToP.y),
//			(argMax == 2) * sgn(cToP.z)
//		);
//	}
//
//
//	float samplePoint(Pcg& gen, Vec3& point_ret) {
//		return 0.0f;
//	}
//};


class Triangle : public Obj {
public:
	Vec3 a;
	Vec3 b;
	Vec3 c;
	Vec3 n;

	float rcp_area;

	Material material;

	Triangle(Vec3 a, Vec3 b, Vec3 c) : a(a), b(b), c(c), n((b-a).cross(c-b).normalized()), rcp_area(2.0f / (b - a).cross(c - b).norm()) {
	}

	float intersect(Ray& ray) {
		Vec3 ba = b - a;
		Vec3 ca = c - a;
		Vec3 roa = ray.o - a;
		Vec3 n = ba.cross(ca);
		Vec3 q = roa.cross(ray.d);
		float d = 1.0f / ray.d.dot(n);
		float u = d * -q.dot(ca);
		float v = d * q.dot(ba);
		if (u<0.0f || u>1.0f || v<0.0f || (u + v)>1.0f) return INFINITY;
		return d * -n.dot(roa);
	}

	Vec3 normal(Vec3 p) {
		return n;
	}

	float samplePoint(Pcg& gen, Vec3& point_ret) const {
		float u[2] = { gen.Uniform(), gen.Uniform() };

		float su0 = std::sqrt(u[0]);
		float b0 = 1.0f - su0;
		float b1 = u[1] * su0;
		point_ret = b0 * a + b1 * b + (1.0f - b0 - b1) * c;

		return rcp_area;
	}

	float lowDiscrepancySamplePoint(Pcg& gen, Vec3& point_ret) {
		// https://pharr.org/matt/blog/2019/03/13/triangle-sampling-1.5
		uint32_t uf = gen.Uniform() * (1ull << 32);  // Fixed point
		float cx = 0.0f;
		float cy = 0.0f;
		float w = 0.5f;

		for (int i = 0; i < 16; i++) {
			uint32_t uu = uf >> 30;
			bool flip = (uu & 3) == 0;

			cy += ((uu & 1) == 0) * w;
			cx += ((uu & 2) == 0) * w;

			w *= flip ? -0.5f : 0.5f;
			uf <<= 2;
		}

		float b0 = cx + w / 3.0f, b1 = cy + w / 3.0f;
		point_ret = b0*a + b1*b +  (1 - b0 - b1)*c;

		return rcp_area;
	}
};

class Mesh : public Obj {
public:
	std::vector<Triangle> triangles;

	// Per triangle material index
	std::vector<unsigned> material_idxs;
	
	std::vector<Material> materials;

	Mesh() {
	}

	void add_triangle(Triangle t, unsigned material_id) {
		triangles.push_back(t);
		material_idxs.push_back(material_id);
	}
	
	void add_material(Material material) {
		materials.push_back(material);
	}

	float intersect(Ray &ray) {
		float min_t = INFINITY;
		Triangle min_triangle = triangles.front();
		for (Triangle triangle : triangles) {
			float t = triangle.intersect(ray);
			if (EPS < t && t < min_t) {
				min_t = t;
				min_triangle = triangle;
			}
		}
		return min_t;
	}

	float samplePoint(Pcg& gen, Vec3& point_ret) const {
		return 0.0f;
	}
};

class EmbreeMesh : public Obj {
public:
	//TODO replace with spans
	Vec3* vertices;
	unsigned* triangles;
	RTCGeometry rtcmesh;

	// Per triangle material index
	std::vector<unsigned> material_idxs;

	std::vector<Material> materials;

	EmbreeMesh(RTCDevice rtcdevice, RTCScene rtcscene, Mesh *mesh) {
		rtcmesh = rtcNewGeometry(rtcdevice, RTC_GEOMETRY_TYPE_TRIANGLE);
		vertices = (Vec3*)rtcSetNewGeometryBuffer(rtcmesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vec3), 3 * mesh->triangles.size());
		triangles = (unsigned*)rtcSetNewGeometryBuffer(rtcmesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), mesh->triangles.size());
		for (unsigned i = 0; i < mesh->triangles.size(); i++) {
			vertices[3 * i] = mesh->triangles[i].a;
			vertices[3 * i + 1] = mesh->triangles[i].b;
			vertices[3 * i + 2] = mesh->triangles[i].c;

			triangles[3 * i] = 3 * i;
			triangles[3 * i + 1] = 3 * i + 1;
			triangles[3 * i + 2] = 3 * i + 2;

		}
		material_idxs = mesh->material_idxs;

		rtcCommitGeometry(rtcmesh);
		unsigned int geomID = rtcAttachGeometry(rtcscene, rtcmesh);
		rtcReleaseGeometry(rtcmesh);

		materials = mesh->materials;
	}

	float intersect(Ray& ray) {
		return 0.0f;
	}

	Vec3 normal(Vec3 p) {
		return Vec3();
	}
	float samplePoint(Pcg& gen, Vec3& point_ret) const {
		return 0.0f;
	}
};

class Sphere : public Obj {
public:
	Vec3 c;
	float r;
	Material material;

	Sphere(Vec3 c, float r) : c(c), r(r) {
	}
	float intersect(Ray& ray) {
		float B = 2.0f * ray.d.dot(ray.o - c);
		float C = (ray.o - c).sqrNorm() - r * r;

		float sqrDiscr = B * B - 4.0f * C;
		if (sqrDiscr < 0.0f) {
			return INFINITY;
		}
		float discr = std::sqrt(sqrDiscr);
		float t1 = (-B - discr) / 2.0f;
		return t1 > EPS ? t1 : (-B + discr) / 2.0f;
	}
	Vec3 normal(Vec3 p) {
		return (p - c) / r;
	}
	float samplePoint(Pcg& gen, Vec3& point_ret) const {
		return 0.0f;
	}
};

//int filter_extra(std::istringstream &s) {
//	std::string x;
//	if (s >> x) {
//		return std::stoi(x.substr(0, x.find_first_of('/')));
//	} else {
//		return -1;
//	}
//
//}

//Mesh *load_mesh(std::string file_name) {
//	std::ifstream ifile(file_name);
//	if (errno) {
//		throw std::system_error(errno, std::system_category(), "Failed to open: " + file_name);
//	} else {
//		std::cout << "Successfully opened file " + file_name << std::endl;
//	}
//	Mesh *mesh = new Mesh();
//
//	std::vector<Vec3> vertices;
//	while (!ifile.eof()) {
//		if (isspace(ifile.peek())) {
//			ifile.get();
//			continue;
//		}
//		std::string line;
//		std::getline(ifile, line);
//		if (line.size() < 2) { // Too short
//		} else if (line.rfind("#", 0) == 0) { // comment
//		} else if (line.rfind("v ", 0) == 0) { // vertex
//			Vec3 vertex;
//			sscanf_s(line.c_str(), "v %f %f %f", &vertex.x, &vertex.y, &vertex.z);
//			vertices.push_back(vertex);
//		} else if (line.rfind("f ", 0) == 0) {//face
//			std::istringstream line_ss(line);
//			std::string word;
//			line_ss >> word;
//
//			int a = filter_extra(line_ss);
//			int b = filter_extra(line_ss);
//			int c = filter_extra(line_ss);
//
//			mesh->add_triangle(Triangle(vertices.at(a - 1), vertices.at(b - 1), vertices.at(c - 1)));
//			//std::cout << mesh->triangles.back().n;
//			b = c;
//			while ((c = filter_extra(line_ss)) >= 0) {
//				mesh->add_triangle(Triangle(vertices.at(a - 1), vertices.at(b - 1), vertices.at(c - 1)));
//				//std::cout << " " << mesh->triangles.back().n;
//				b = c;
//			}
//			//std::cout << std::endl;
//		}
//	}
//	return mesh;
//}
