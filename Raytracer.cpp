// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <thread>
#include <functional>
#include <string.h>
#include <array>

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <random> // TODO Replace with QMC generators

#include <chrono>
#include <ctime>  

#include <omp.h>
#include <SFML/Graphics.hpp>
#include <embree3/rtcore.h>

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"

#define WIDTH 800
#define HEIGHT 800
#define FOVX 60.0f
#define ACTIVE_SCENE 3

#define AMBIENT 0.3f
#define SAMPLES_PER_PIXEL 4
#define THREADS 2

// How much emission a material should have to be sampled by NEE
#define NEE_EMISSION_THRESHOLD 0.05

// Unused
#define MIN_BOUNCES 2
#define MAX_BOUNCE_PROB 0.98f
#define HEMISPHERE_AREA (M_PI*2.0f)

using namespace std;

class Img {
public:
	array<array<Vec3, WIDTH>, HEIGHT> *pixels;
	unsigned int sample_count;
	Img() {
		pixels = new array<array<Vec3, WIDTH>, HEIGHT>();
		clear();
	}
	
	void clear() {
		for (int py = 0; py < HEIGHT; py++) {
			for (int px = 0; px < WIDTH; px++) {
				(*pixels)[py][px] = Vec3();
			}
		}
		sample_count = 0;
	}

	void update(int px, int py, Vec3 color) {
		// TODO use Kahan summation
		(*pixels)[py][px] = ((*pixels)[py][px] * (float) (sample_count) + color) / ((float) (sample_count + 1));
	}
};

struct Intersection {
	float t;
	Obj * hitObj;
	Material material;
	Vec3 normal;
};

void rtcErrorFunction(void* userPtr, enum RTCError error, const char* str) {
	printf("RTCERROR %d: %s\n", error, str);
}

class Scene {
public:
	vector<Obj *> objs;
	vector<Mesh *> meshes;

	vector<Triangle> lights;

	vector<Material> materials;
	
	RTCDevice rtcdevice;
	RTCScene rtcscene;

	Scene() {
		rtcdevice = rtcNewDevice(nullptr);
		rtcSetDeviceErrorFunction(rtcdevice, rtcErrorFunction, NULL);
		rtcscene = rtcNewScene(rtcdevice);
		rtcSetSceneBuildQuality(rtcscene, RTC_BUILD_QUALITY_HIGH);
	}

	~Scene() {
		rtcReleaseScene(rtcscene);
		rtcReleaseDevice(rtcdevice);
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

		//if (!reader.Warning().empty()) {
		//	std::cout << "TinyObjReader: " << reader.Warning();
		//}

		Mesh *mesh = new Mesh();

		auto& attrib = reader.GetAttrib();
		auto& shapes = reader.GetShapes();
		for (size_t i = 0; i < reader.GetMaterials().size(); i++) {
			mesh->add_material({ Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) });
		}

		//TODO rewrite to just use a basic buffer
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
					mesh->add_triangle(Triangle(a, b, c), shapes[s].mesh.material_ids[f]);
					b = c;
				}
				index_offset += fv;

			}
		}
		return mesh;
	}

	vector<EmbreeMesh> e_meshes;
	void addMesh(Mesh *mesh) {
		meshes.push_back(mesh);
		e_meshes.push_back(EmbreeMesh(rtcdevice, rtcscene, mesh));
	}

	void commit() {
		rtcCommitScene(rtcscene);

		for (const EmbreeMesh& mesh : e_meshes) {
			for (size_t i = 0; i < mesh.material_idxs.size(); i++) {
				if (mesh.materials[mesh.material_idxs[i]].emission > NEE_EMISSION_THRESHOLD) {
					lights.push_back(
						Triangle(
							mesh.vertices[3 * mesh.triangles[i]],
							mesh.vertices[3 * mesh.triangles[i] + 1],
							mesh.vertices[3 * mesh.triangles[i] + 2]
					));
				}
			}
		}
	}

	void castRay(Ray& ray, Intersection& intersection) {
		RTCIntersectContext rtccontext;
		rtcInitIntersectContext(&rtccontext);
		RTCRayHit rayhit;
		rayhit.ray.org_x = ray.o.x;
		rayhit.ray.org_y = ray.o.y;
		rayhit.ray.org_z = ray.o.z;

		rayhit.ray.dir_x = ray.d.x;
		rayhit.ray.dir_y = ray.d.y;
		rayhit.ray.dir_z = ray.d.z;

		rayhit.ray.tnear = EPS;
		rayhit.ray.tfar = std::numeric_limits<float>::infinity();

		rayhit.ray.flags = 0;
		rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

		rtcIntersect1(rtcscene, &rtccontext, &rayhit);
		// Use rayhit.hit.geomid to determine the hit object, primId to determine the face
		if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
			auto mesh = meshes[rayhit.hit.geomID];
			intersection.hitObj = &(mesh->triangles[rayhit.hit.primID]);
			intersection.material = mesh->materials[mesh->material_idxs[rayhit.hit.primID]];

			intersection.normal = Vec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z).normalized();
		}
		intersection.t = rayhit.ray.tfar;
	}

	bool isOccluded(Ray& ray, float distance) {
		if (distance < 2.0f*EPS) {
			return false;
		}

		RTCIntersectContext rtccontext;
		rtcInitIntersectContext(&rtccontext);

		RTCRay rtcray;
		rtcray.org_x = ray.o.x;
		rtcray.org_y = ray.o.y;
		rtcray.org_z = ray.o.z;

		rtcray.dir_x = ray.d.x;
		rtcray.dir_y = ray.d.y;
		rtcray.dir_z = ray.d.z;

		rtcray.tnear = EPS;
		rtcray.tfar = distance - EPS;

		rtcray.flags = 0;

		rtcOccluded1(rtcscene, &rtccontext, &rtcray);

		return rtcray.tfar == -std::numeric_limits<float>::infinity();
	}

	float sampleRandomPointOnLight(mt19937& gen, Vec3 &point_ret) {
		uniform_int_distribution<unsigned> uniformLight(0, lights.size()-1);
		float pdf = lights[uniformLight(gen)].samplePoint(gen, point_ret);
		return pdf;
	}
};

Scene scene;

uniform_real_distribution<float> pix(-0.5f, 0.5f);
uniform_real_distribution<float> uniformAngle(0.0f, 2.0f*M_PI);

Obj *light = nullptr;
Material redMat = { Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) };
void buildScene(int i) {
	Material mirrorMat = { Vec3(1.0f), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(diffuse) };
	//Material redMat = { Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) };
	Material greenMat = { Vec3(0.12f, 0.45f, 0.15f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };

	// Using specular light sauces creates a lot of noise
	Material lightMat = { Vec3(1.0f), 16.0f, Surface(specular) };
	Material ovenMat = { Vec3(0.5f), 0.5f, Surface(diffuse) };

	scene.addMesh(scene.load_mesh("geometry/CornellBox-Original.obj"));
	scene.meshes[0]->materials[7].emission = 4.0f;
	scene.meshes[0]->materials[6] = specularMat;

	//scene.meshes.push_back(load_mesh("geometry/teapot.obj"));
	//scene.meshes.front()->material = diffuseMat;
	//switch (i) {
	//default:
	//	break;
	//case 0: // Plane
	//	scene.addObject(new Plane(Vec3(0.0f, 2.05f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
	//	//scene.addObject(new Sphere(Vec3(2.0f, 1.0f, -4.0f), 0.8f), lightMat);
	//	scene.addObject(new Sphere(Vec3(1.7f, 0.5f, -4.0f), 1.3f), specularMat);
	//	scene.addObject(new Box(Vec3(-1.4f, 0.5f, -3.0f), Vec3(-0.4f, 1.5f, -2.0f)), specularMat);
	//	break;
	//case 1: // Cornell box
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), greenMat);
	//	scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), redMat);

	//	scene.addObject(new Sphere(Vec3(400.0f, 440.0f, -600.0f), 200.3f), specularMat);
	//	scene.addObject(new Box(Vec3(-200.0f, -800.0f, -500.0f), Vec3(200.0f, -750.0f, -100.0f)), specularMat);

	//	scene.addObject(new Box(Vec3(-200.0f, -800.0f, -500.0f), Vec3(200.0f, /*-795.0f*/ -799.9f, -100.0f)), lightMat);
	//	light = scene.objs.back();
	//	break;
	//case 2: // Oven test
	//	/*
	//	* The oven test is any encosed room with surface emission 0.5 and reflectance 0.5. So we expect a pixel value
	//	* 0.5*(0.5 + 0.5(0.5 + 0.5(...)) = 1.
	//	*/
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), ovenMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), ovenMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), ovenMat);
	//	scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), ovenMat);
	//	scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), ovenMat);
	//	scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), ovenMat);
	//	break;
	//case 3:
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), diffuseMat);
	//	scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), greenMat);
	//	scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), redMat);

	//	scene.addObject(new Sphere(Vec3(0.0f, -040.0f, -500.0f), 200.3f), diffuseMat);
	//	scene.addObject(new Box(Vec3(-50.0f, -800.0f, -650.0f), Vec3(50.0f, /*-795.0f*/ -799.9f, -550.0f)), lightMat);
	//	light = scene.objs.back();
	//}
	//scene.addObject(new Box(Vec3(1.0f, 1.0f, -3.0f), Vec3(1.0f, 1.0f, -3.0f) + Vec3(1.0f, 1.0f, 0.0f)), diffuseMat);
	/*

	/*for (Obj *obj : scene) {
		if (obj->emission > 0.0f) {
			lights.push_back(obj);
		}
	}*/

	scene.commit();
}

// Transform pixel coordinates to perspective rays
// Camera is at (0,0,0) facing (0,0,-1)
void cameraRay(float px, float py, Ray& ray) {
	float x = (2.0f * px - WIDTH) / WIDTH * tanf(FOVX * M_PI / 180.0f / 2.0f);
	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	float y = (2.0f * py - HEIGHT) / HEIGHT * tanf(((float) HEIGHT) / WIDTH * FOVX * M_PI / 180.0f / 2.0f);
	float z = -1.0f;
	ray.o = Vec3(0.0f,1.0f,3.0f);
	ray.d = Vec3(x,-y,z).normalized();

	//return Ray(Vec3(), Vec3(x, y, z));
}

Vec3 cosineSampleHemisphere(mt19937 &gen) {
	float r_sqr = uniform01(gen);
	float phi = uniformAngle(gen);
	float r = sqrt(r_sqr);
	float x = cos(phi)*r;
	float y = sin(phi)*r;
	return Vec3(x, y, sqrt(1.0f - x*x - y*y));
}
Vec3 uniformSampleHemisphere(mt19937 &gen) {
	float z = uniform01(gen);
	float theta = uniformAngle(gen);
	float r = sqrt(1.0f-z*z);
	return Vec3(cos(theta)*r, sin(theta)*r, z);
}

// Creates an orthonormal basis
void onb(Vec3 v1, Vec3 &v2, Vec3 &v3) {
	if (abs(v1.x) > abs(v1.y)) {
		v2 = Vec3(-v1.z, 0.0f, v1.x) / sqrt(v1.z*v1.z + v1.x*v1.x);
	} else {
		v2 = Vec3(0.0f, v1.z, -v1.y) / sqrt(v1.z*v1.z + v1.y*v1.y);
	}
	v3 = v1.cross(v2);
}

// cos_t is the dot product of the normal and the incident vectors
float schlickApprox(float r, float cos_t) {
	float R0 = (r - 1.0f) / (r + 1.0f);
	R0 *= R0;
	float x = 1.0f - cos_t;
	float x2 = x * x;
	return R0 + (1.0f - R0)*x2*x2*x;
}

Vec3 henyey_greenstein(Vec3 in, float g, mt19937 &gen) {
	float sqr_part = (1.0f - g * g) / (1.0f + g * g + 2.0f*g*uniform01(gen));
	float cos_t = -(1.0f + g * g - (sqr_part*sqr_part)) / 2.0f / g;
	float sin_t = sqrtf(max(0.0f, 1.0f - cos_t*cos_t));
	float phi = uniformAngle(gen);
	Vec3 v2, v3;
	onb(in, v2, v3);
	//v3.normalized();
	//cout << in.norm() << endl;
	return v2 * sin_t*cos(phi) + v3 * sin_t*sin(phi) + in * cos_t;
}
float henyey_greenstein_p(float cos_t, float g=1.0f) {
	return 1.0f/4.0f/M_PI*(1.0f-g*g)/(1.0f+g*g+2.0f*g*powf(cos_t,1.5f));
}

Vec3 direction(float theta, float phi) {
	return Vec3(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
}

int c = 0;
Vec3 rayTrace(Ray& ray, mt19937 &gen) {
	Vec3 attenuation = Vec3(1.0f);
	Vec3 color = Vec3(0.0f);
	Obj *last = nullptr;
	//int i = 0;
	//int a = 0;
	float mis_brdf_pdf = -1.0f; // Negative when mis was not used
	while (true) {
		Intersection intersection;
		scene.castRay(ray, intersection);
		
		//TODO: Weigh this by MIS
		if (isinf(intersection.t)) {
			float a = ray.d.dot(Vec3(0.2f, -0.8f, -0.4f).normalized());
			/*if (a > 0.999) {

				color = color + attenuation*Vec3(5);
			} else if (a > 0.96) {
				color = color + attenuation*Vec3(5)*(a-0.96f)*(a - 0.96f)/ (0.999f - 0.96f) /(0.999f-0.96f);
			}*/
			//color = color + attenuation*Vec3(0.5f,0.70f,0.8f);
			break;
		}
		float density = 0.0f;

		Vec3 hp = intersection.t*ray.d + ray.o;

		Obj *obj = intersection.hitObj;
		Vec3 n = intersection.normal;
		
		Material material = intersection.material;
		//obj->material = redMat;
		// TODO: remove this
		//obj->material.albedo = n*0.7f + 0.3f;
		//obj->material.albedo.x = min(abs(obj->material.albedo.x), 1.0f);
		//obj->material.albedo.y = min(abs(obj->material.albedo.y), 1.0f);
		//obj->material.albedo.z = min(abs(obj->material.albedo.z), 1.0f);

		Vec3 emission = Vec3(material.emission);

		if (emission.max() > 0.0f) {
			// Apply MIS
			if (mis_brdf_pdf > 0.0f) {
				float solidAngle = abs(-ray.d.dot(n)) * 100 * 100 / intersection.t / intersection.t;
				// We have 2 multiplications of mis_brdf_pdf because, we have already divided the attenuation by mis_brdf_pdf
				color = color + emission*attenuation *
					mis_brdf_pdf * mis_brdf_pdf / (1 / solidAngle / solidAngle + mis_brdf_pdf * mis_brdf_pdf);
					//* exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER
				
			} else {
				color = color + emission * attenuation;
					//* exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER
			}
		}

		// Medium intersection
		//float medium_t = -log(1.0f - uniform01(gen)) / density;
		float medium_t = INFINITY;
		if (intersection.t > medium_t) {
			//hit medium

			ray.o = medium_t * ray.d + ray.o;

			// MIS direct light
			Vec3 sampledLight;
			float pdf = light->samplePoint(gen, sampledLight);
			ray.d = (sampledLight - ray.o).normalized();
			Intersection v;

			// TODO change to use occluded check, technically light == v.hitobj is incorrect
			scene.castRay(ray, v);
			if (light == v.hitObj) {
				float pl = 1.0f / 100.0f / 100.0f;
				float solidAngle = abs(-ray.d.dot(v.normal)) / pl / v.t / v.t;

				color = color + v.material.emission * attenuation * material.albedo / M_PI * abs(n.dot(ray.d)) *
					1 / solidAngle / (1 / solidAngle / solidAngle + abs(ray.d.dot(n)) / M_PI * abs(ray.d.dot(n)) / M_PI) *
					exp(-density * v.t); // MEDIUM TERM REMOVE LATER
			}

			Vec3 d = ray.d;
			ray.d = henyey_greenstein(ray.d, 0.1f, gen);
			attenuation = attenuation * Vec3(0.999);
			mis_brdf_pdf = henyey_greenstein_p(ray.d.dot(d), 0.1f);
			continue;
		}
		ray.o = hp;

		// TODO
		// Some rays get eliminated before they can even do anything,
		// So either try setting minimum bounces or rearranging some of the code
		float bounce_probability = min(attenuation.max(), MAX_BOUNCE_PROB);
		
		if (uniform01(gen) > bounce_probability) {
			break;
		}
		attenuation = attenuation / bounce_probability;

		if (material.surface == reflective) {
			float cos_t = ray.d.dot(n);

			ray.d = (ray.d - (n * cos_t * 2.0f));
			attenuation = attenuation * material.albedo;
			mis_brdf_pdf = -1.0f;
		} else if (material.surface == diffuse) {
			// TODO actually do this properly
			float pl = 1.0f / 100.0f / 100.0f;
			//if (light != obj) {
			//	Vec3 sampledLight = light->samplePoint(gen);
			//	ray.d = (sampledLight - ray.o).normalized();
			//	Intersection v = scene.castRay(ray);
			//	if ( light == v.hitObj) {
			//		float solidAngle = abs(-ray.d.dot(v.hitObj->normal(ray.o + v.t*ray.d)))/pl / v.t / v.t;

			//		color = color + light->material.emission * attenuation * material.albedo / M_PI * abs(n.dot(ray.d)) *
			//			1 / solidAngle /( 1/solidAngle/solidAngle + abs(ray.d.dot(n))/M_PI* abs(ray.d.dot(n)) / M_PI) *
			//			exp(-density*v.t); // MEDIUM TERM REMOVE LATER
			//	}
			//}
			//color = color + emission * attenuation;

			// TODO: Use onb's instead
			Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0f, 0.0f, 1.0f));
			ray.d = rotMatrix * cosineSampleHemisphere(gen);
			float cos_t = abs(ray.d.dot(n));
			//mis_brdf_pdf = cos_t / M_PI;
			attenuation = attenuation * material.albedo;/**cos_t/(cos_t/M_PI )*/ // * pi (Surface area) / (pi (lambertian albedo constant))
		} else if (material.surface == specular) {
			// TODO Be able to handle materials with different refractive indexes
			float r = 1.0f / 1.5f;
			float cos_t1 = -n.dot(ray.d);
			if (cos_t1 < 0.0f) {
				// We're inside the specular object
				cos_t1 *= -1;
				r = 1.0f / r;
				n = -n;
			}

			// Check for total internal reflection, then choose refraction based on fresnel

			//TODO: See how much slower computing the real fresnel reflectance is. (will be twice as likely to call sqrt() than usual, but will be more accurate)
			if (1.0f - r * r*(1.0f - cos_t1 * cos_t1) > 0.0f &&
				uniform01(gen) > ((r > 1.0f) ? schlickApprox(r, cos_t1) : schlickApprox(r, sqrt(1.0f - r * r*(1.0f - cos_t1 * cos_t1))))) {
				// Refraction through the specular surface
				float cos_t2 = sqrt(1.0f - r * r*(1.0f - cos_t1 * cos_t1));
				ray.d = (r*ray.d + (r*cos_t1 - cos_t2)*n);
			} else {
				// Reflection off the specular surface
				ray.d = (ray.d + (n * cos_t1 * 2));
			}
			attenuation = attenuation * material.albedo;
			mis_brdf_pdf = -1.0f;
		}
		// Some diagnostic tools
		// i++;
		//c += last == obj && obj == light;
		//a += last == obj && obj == light;
		//last = obj;
	}
	//if (i > 1) cout << color << " " << i << endl;
	//if (a > 1) cout << color << " " << a << endl;
	//if (color.max() < 0.4) cout << color << endl;
	return color;
}

random_device rd;
mt19937 gens[THREADS];

void render(Img &img) {
	auto t1 = chrono::high_resolution_clock::now();
	if (img.sample_count == 0) {
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int py = 0; py < HEIGHT; py++) {
			mt19937& gen = gens[omp_get_thread_num()];
			for (int px = 0; px < WIDTH; px++) {
				Ray ray = Ray(Vec3(), Vec3());
				Vec3 colorVec = Vec3();
				for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
					cameraRay(px + pix(gen), py + pix(gen), ray);
					colorVec = colorVec + rayTrace(ray, gen);
				}
				(*img.pixels)[py][px] = colorVec/SAMPLES_PER_PIXEL;
			}
		}
	} else {
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int py = 0; py < HEIGHT; py++) {
			mt19937& gen = gens[omp_get_thread_num()];
			for (int px = 0; px < WIDTH; px++) {
				Ray ray = Ray(Vec3(), Vec3());
				Vec3 colorVec = Vec3();
				for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
					cameraRay(px + pix(gen), py + pix(gen), ray);
					colorVec = colorVec + rayTrace(ray, gen);
				}
				img.update(px, py, colorVec/SAMPLES_PER_PIXEL);
			}
		}

	}
		
	// Performance metrics
	float seconds = (chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0);
	printf("Took %f seconds\n", seconds);
	printf("Cast %d rays\n", SAMPLES_PER_PIXEL * WIDTH*HEIGHT);
	printf("%f samples per second \n", SAMPLES_PER_PIXEL * WIDTH*HEIGHT / seconds);
}

void toneMap(array<array<Vec3, WIDTH>, HEIGHT> *img) {
	//for (int py = 0; py < HEIGHT; py++) {
	//	for (int px = 0; px < WIDTH; px++) {
	//		(*img)[py][px] = Vec3((*img)[py][px].x, (*img)[py][px].y, (*img)[py][px].z);
	//	}
	//}
}


string getDateTime() {
	// get the current time
	time_t t = time(0);

	tm now;
	localtime_s(&now, &t);
	char buffer[15];
	std::strftime(buffer, 32, "%Y%m%d%H%M%S", &now);
	return string(buffer);
}

void render_loop(Img &img) {
	auto t0 = chrono::high_resolution_clock::now();
	for (int i = 0; i < THREADS; i++) {
		gens[i] = mt19937(hash<int>{}(i));
	}
	while (true) {
		render(img);
		img.sample_count++;

		float seconds = (chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t0).count() / 1000000.0);
		printf("Accumulated %d samples per pixel over %f seconds\n\n", SAMPLES_PER_PIXEL*img.sample_count, seconds);
	}
}

void process_image(sf::Image &image, array<array<Vec3, WIDTH>, HEIGHT> *img) {
	toneMap(img);
	// Load image
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			image.setPixel(px, py, (*img)[py][px].toColor());
		}
	}
}

void gui_thread(Img &img) {
	
	// Tone map resulting image

	sf::Image image;
	image.create(WIDTH, HEIGHT);

	//process_image(image, img);

	sf::ContextSettings contextSettings;
	contextSettings.sRgbCapable = true;

	// Init texture
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer", sf::Style::Default, contextSettings);
	sf::Texture texture;
	sf::Sprite sprite;
	texture.setSrgb(false);

	ofstream ofile;
	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			this_thread::sleep_for(100ms);
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				switch (event.key.code) {
				case sf::Keyboard::S:
					printf("Saving render...\n");
					image.saveToFile(".\\renders\\render" + getDateTime() + ".png"); // TODO save the file as srgb
					ofile.open(".\\renders\\renderData" + getDateTime() + ".txt");
					for (int py = 0; py < HEIGHT; py++) {
						for (int px = 0; px < WIDTH; px++) {
							Vec3 v = (*img.pixels)[py][px];
							ofile << (int)v.x << ',' << (int)v.y << ',' << (int)v.z << ' ';
						}
						ofile << endl;
					}
					ofile.close();

					printf("Render saved\n");
					break;
				case sf::Keyboard::Escape:
					window.close();
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		}
		window.clear();
		window.draw(sprite);
		window.display();

		process_image(image, img.pixels);
		texture.loadFromImage(image);
		sprite.setTexture(texture);
	}
}

int main() {
	if (!std::numeric_limits<float>::is_iec559) {
		printf("Machine architecture must implement IEEE 754.\n");
		return 0;
	}
	
	auto t1 = chrono::high_resolution_clock::now();
	buildScene(ACTIVE_SCENE);
	float seconds = (chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0);
	printf("Scene built in %f seconds\n", seconds);

	Img img = Img();
	img.clear();

	// Render scene
	thread render_thread = thread(render_loop, ref(img));
	gui_thread(img);
	exit(0);
	return 0;
}