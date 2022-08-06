#pragma once

#include <embree3/rtcore.h>

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"
#include "Camera.h"
#include "Scene.h"

struct Intersection {
	float t;
	Obj* hitObj;
	Material material;
	Vec3 normal;
};

void rtcErrorFunction(void* userPtr, enum RTCError error, const char* str) {
	printf("RTCERROR %d: %s\n", error, str);
}

class Scene {
public:
	std::vector<Obj*> objs;
	std::vector<Mesh*> meshes;

	std::vector<Triangle> lights;

	std::vector<Material> materials;

	std::vector<EmbreeMesh> e_meshes;

#ifdef SPHERES
	std::vector<Sphere> spheres;
#endif

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

	Mesh* load_mesh(std::string file_name) {

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

		Mesh* mesh = new Mesh();

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

				// Tesselate the face
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

	void addMesh(Mesh* mesh) {
		meshes.push_back(mesh);
		e_meshes.push_back(EmbreeMesh(rtcdevice, rtcscene, mesh));
	}

	void commit() {
		rtcCommitScene(rtcscene);
		for (size_t i = 0; i < meshes.size(); i++)
		{
			e_meshes[i].materials = meshes[i]->materials;
		}

		for (const EmbreeMesh& mesh : e_meshes) {
			for (size_t i = 0; i < mesh.material_idxs.size(); i++) {
				if (mesh.materials[mesh.material_idxs[i]].emission > NEE_EMISSION_THRESHOLD) {
					lights.push_back(
						Triangle(
							mesh.vertices[3 * mesh.triangles[i]],
							mesh.vertices[3 * mesh.triangles[i] + 1],
							mesh.vertices[3 * mesh.triangles[i] + 2]
						));
					lights.back().material = mesh.materials[mesh.material_idxs[i]];
				}
			}
		}
	}

	void castRay(Ray& ray, Intersection& intersection) const {
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
#ifdef SPHERES
		intersection.t = std::numeric_limits<float>::infinity();
		for (Sphere sphere : spheres) {
			float t = sphere.intersect(ray);
			if (t > EPS && t < intersection.t) {
				intersection.t = t;
				intersection.hitObj = &sphere;
				intersection.material = sphere.material;
				intersection.normal = sphere.normal(ray.o + t * ray.d);
			}
		}

		if (intersection.t < rayhit.ray.tfar) {
			return;
		}
#endif

		// Use rayhit.hit.geomid to determine the hit object, primId to determine the face
		if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
			auto mesh = meshes[rayhit.hit.geomID];
			intersection.hitObj = &(mesh->triangles[rayhit.hit.primID]);
			intersection.material = mesh->materials[mesh->material_idxs[rayhit.hit.primID]];

			intersection.normal = Vec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z).normalized();
		}
		intersection.t = rayhit.ray.tfar;
	}

	bool isOccluded(Ray& ray, float distance) const {
		if (distance < 2.0f * EPS) {
			return false;
		}
#ifdef SPHERES
		for (Sphere sphere : spheres) {
			float t = sphere.intersect(ray);
			if (t > EPS && t < distance - EPS) {
				return true;
			}
		}
#endif

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

		// Negative infinity indicates a hit
		return rtcray.tfar == -std::numeric_limits<float>::infinity();
	}

	float sampleRandomPointOnLight(Pcg& gen, Vec3& point_ret, Triangle const** light) const {
		if (lights.size() == 0) {
			return 0.0f;
		}
		// TODO implement integer sampling, I don't think this is even threadsafe
		*light = &(lights[rand() % 2]);
		float pdf = (*light)->samplePoint(gen, point_ret);
		return pdf / lights.size();
	}
};