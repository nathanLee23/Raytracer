#pragma once

#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

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

	Triangle(Vec3 a, Vec3 b, Vec3 c) : a(a), b(b), c(c), n((b-a).cross(b-c).normalized()) {
	}

	float intersect(Ray& ray) {
		Vec3 ba = b - a;
		Vec3 ca = c - a;
		Vec3 roa = ray.o - a;
		Vec3 n = ba.cross(ca);
		Vec3 q = roa.cross(ray.d);
		float d = 1.0 / ray.d.dot(n);
		float u = d * -q.dot(ca);
		float v = d * q.dot(ba);
		if (u<0.0 || u>1.0 || v<0.0 || (u + v)>1.0) return INFINITY;
		return d * -n.dot(roa);
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
		bounds.max.x = std::max({ bounds.max.x, t.a.x, t.b.x, t.c.x });
		bounds.max.y = std::max({ bounds.max.y, t.a.y, t.b.y, t.c.y });
		bounds.max.z = std::max({ bounds.max.z, t.a.z, t.b.z, t.c.z });

		bounds.min.x = std::min({ bounds.min.x, t.a.x, t.b.x, t.c.x });
		bounds.min.y = std::min({ bounds.min.y, t.a.y, t.b.y, t.c.y });
		bounds.min.z = std::min({ bounds.min.z, t.a.z, t.b.z, t.c.z });

	}
	float intersect(Ray &ray) {
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

	float tri_intersect(Ray &ray, Triangle **min_triangle) {
		float min_t = bounds.intersect(ray);
		if (std::isinf(min_t))
			return min_t;
		min_t = INFINITY;
		for (int i = 0; i < triangles.size(); i ++) {
			Triangle triangle = triangles[i];
			float t = triangle.intersect(ray);
			if (t < min_t) {
				min_t = t;
				*min_triangle = &triangles[i];
				// TODO: do something more efficient here
				// like make the material part of the return,
				// and get rid of materials from triangles
				(*min_triangle)->material = material;
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

Mesh *load_mesh(std::string file_name) {
	
	tinyobj::ObjReaderConfig reader_config;
	//reader_config.mtl_search_path = "./"; // Path to material files

	tinyobj::ObjReader reader;

	if (!reader.ParseFromFile(file_name, reader_config)) {
		if (!reader.Error().empty()) {
			std::cerr << "TinyObjReader: " << reader.Error();
		}
		exit(1);
	}

	if (!reader.Warning().empty()) {
		std::cout << "TinyObjReader: " << reader.Warning();
	}

	auto& attrib = reader.GetAttrib();
	auto& shapes = reader.GetShapes();
	Mesh *mesh = new Mesh();
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			int fv = shapes[s].mesh.num_face_vertices[f];
			
			tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + 0];
			tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
			tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
			tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
			Vec3 a = Vec3(vx, vy, vz);

			idx = shapes[s].mesh.indices[index_offset + 1];
			vx = attrib.vertices[3 * idx.vertex_index + 0];
			vy = attrib.vertices[3 * idx.vertex_index + 1];
			vz = attrib.vertices[3 * idx.vertex_index + 2];
			Vec3 b = Vec3(vx, vy, vz);

			// Loop over vertices in the face.
			for (size_t v = 2; v < fv; v++) {
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
				tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
				tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
				Vec3 c = Vec3(vx, vy, vz);
				mesh->add_triangle(Triangle(a,b,c));
				b = c;
			}
			index_offset += fv;

			// per-face material
			//shapes[s].mesh.material_ids[f];
		}
	}
	return mesh;
}

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
