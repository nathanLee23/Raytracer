#pragma once

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"
#include "Camera.h"
#include "Scene.h"

#define MAX_BOUNCE_PROB 0.99f

// Unused
#define MIN_BOUNCES 2
#define HEMISPHERE_AREA (M_PI*2.0f)

Vec3 cosineSampleHemisphere(Pcg& gen) {
	float r_sqr = gen.Uniform();
	float phi = gen.Uniform() * 2.0f * M_PI;
	float r = sqrt(r_sqr);
	float x = cos(phi) * r;
	float y = sin(phi) * r;
	return Vec3(x, y, sqrt(1.0f - x * x - y * y));
}

Vec3 uniformSampleHemisphere(Pcg& gen) {
	float z = gen.Uniform();
	float theta = gen.Uniform() * 2.0f * M_PI;
	float r = sqrt(1.0f - z * z);
	return Vec3(cos(theta) * r, sin(theta) * r, z);
}

// Creates an orthonormal basis, assumes v1 has unit norm
// From T. Duf et al. Building an Orthonormal Basis, Revisited. Journal of Computer Graphics Techniques
void onb(const Vec3 v1, Vec3& v2, Vec3& v3) {
	float sign = copysignf(1.0f, v1.z);
	const float a = -1.0f / (sign + v1.z);
	const float b = v1.x * v1.y * a;
	v2 = Vec3(1.0f + sign * v1.x * v1.x * a, sign * b, -sign * v1.x);
	v3 = Vec3(b, sign + v1.y * v1.y * a, -v1.y);
}

// cos_t is the dot product of the normal and the incident vectors
float schlickApprox(float r, float cos_t) {
	float R0 = (r - 1.0f) / (r + 1.0f);
	R0 *= R0;
	float x = 1.0f - cos_t;
	float x2 = x * x;
	return R0 + (1.0f - R0) * x2 * x2 * x;
}

Vec3 henyey_greenstein(Vec3 in, float g, Pcg& gen) {
	float sqr_part = (1.0f - g * g) / (1.0f + g * g + 2.0f * g * gen.Uniform());
	float cos_t = -(1.0f + g * g - (sqr_part * sqr_part)) / 2.0f / g;
	float sin_t = sqrtf(std::max(0.0f, 1.0f - cos_t * cos_t));
	float phi = gen.Uniform() * 2.0f * M_PI;
	Vec3 v2, v3;
	onb(in, v2, v3);
	//v3.normalized();
	//cout << in.norm() << endl;
	return v2 * sin_t * cos(phi) + v3 * sin_t * sin(phi) + in * cos_t;
}
float henyey_greenstein_p(float cos_t, float g = 1.0f) {
	return 1.0f / 4.0f / M_PI * (1.0f - g * g) / (1.0f + g * g + 2.0f * g * powf(cos_t, 1.5f));
}

Vec3 direction(float theta, float phi) {
	return Vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
}

// The light emitted by the sky from direction
Vec3 skyEmission(Vec3 const direction) {
	float a = direction.dot(Vec3(0.2f, -0.8f, -0.4f).normalized());
	Vec3 color;
	if (a > 0.999) {
		color = Vec3(5);
	}
	else if (a > 0.96) {
		color = Vec3(5) * (a - 0.96f) * (a - 0.96f) / (0.999f - 0.96f) / (0.999f - 0.96f);
	}
	color += Vec3(0.5f, 0.70f, 0.8f);

	return color;
}

Vec3 rayTraceNormals(const Scene& scene, Ray& ray, Pcg& gen) {
	Intersection intersection;
	scene.castRay(ray, intersection);
	if (isinf(intersection.t)) {
		return skyEmission(ray.d);
	}
	else {
		return intersection.normal * 0.5f + 0.5f;
	}
}

Vec3 pathTrace(const Scene& scene, Ray& ray, Pcg& gen) {
	Vec3 attenuation = Vec3(1.0f);
	Vec3 color = Vec3(0.0f);
	Obj* last = nullptr;
	//int i = 0;
	//int a = 0;
	//static int c = 0;
	float mis_brdf_pdf = -1.0f; // Negative when mis was not used
	float nee_pdf;

	while (true) {
		Intersection intersection;
		scene.castRay(ray, intersection);

		//TODO: Weigh this by MIS
		if (isinf(intersection.t)) {
			//color += attenuation*skyEmission(ray.d);
			break;
		}
		float density = 0.0f;

		Vec3 hp = intersection.t * ray.d + ray.o;

		Obj* obj = intersection.hitObj;
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
				float solidAngle = abs(-ray.d.dot(n)) / nee_pdf / intersection.t / intersection.t;
				// TODO, nee_pdf here is incorrect, fix it (should use the sample pdf of the the object we hit)

				// We have 2 multiplications of mis_brdf_pdf because, we have already divided the attenuation by mis_brdf_pdf
				color += emission * attenuation *
					mis_brdf_pdf * mis_brdf_pdf / (1.0f / solidAngle / solidAngle + mis_brdf_pdf * mis_brdf_pdf);
				//* exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER

			}
			else {
				color += emission * attenuation;
				//* exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER
			}
		}

		// Medium intersection
		//float medium_t = -log(1.0f - gen.Uniform()) / density;
		//float medium_t = INFINITY;
		//if (intersection.t > medium_t) {
		//	//hit medium

		//	ray.o = medium_t * ray.d + ray.o;

		//	// MIS direct light
		//	Vec3 sampledLight;
		//	float pdf = light->samplePoint(gen, sampledLight);
		//	ray.d = (sampledLight - ray.o).normalized();
		//	Intersection v;

		//	// TODO change to use occluded check, technically light == v.hitobj is incorrect
		//	scene.castRay(ray, v);
		//	if (light == v.hitObj) {
		//		float pl = 1.0f / 100.0f / 100.0f;
		//		float solidAngle = abs(-ray.d.dot(v.normal)) / pl / v.t / v.t;

		//		color = color + v.material.emission * attenuation * material.albedo / M_PI * abs(n.dot(ray.d)) *
		//			1 / solidAngle / (1 / solidAngle / solidAngle + abs(ray.d.dot(n)) / M_PI * abs(ray.d.dot(n)) / M_PI) *
		//			exp(-density * v.t); // MEDIUM TERM REMOVE LATER
		//	}

		//	Vec3 d = ray.d;
		//	ray.d = henyey_greenstein(ray.d, 0.1f, gen);
		//	attenuation = attenuation * Vec3(0.999);
		//	mis_brdf_pdf = henyey_greenstein_p(ray.d.dot(d), 0.1f);
		//	continue;
		//}
		ray.o = hp;

		// TODO
		// Some rays get eliminated before they can even do anything,
		// So either try setting minimum bounces or rearranging some of the code
		float bounce_probability = std::min(attenuation.max(), MAX_BOUNCE_PROB);

		if (gen.Uniform() > bounce_probability) {
			break;
		}
		attenuation = attenuation / bounce_probability;

		switch (material.surface)
		{
		case reflective: {
			float cos_t = ray.d.dot(n);

			ray.d = (ray.d - (n * cos_t * 2.0f));
			attenuation = attenuation * material.albedo;
			mis_brdf_pdf = -1.0f;
			break;
		}
		case diffuse: {
			bool yes = false;
			// NEE
			if (scene.lights.size() != 0 && material.emission < NEE_EMISSION_THRESHOLD) {
				Triangle const* sampledLight;
				Vec3 sampledPoint;
				nee_pdf = scene.sampleRandomPointOnLight(gen, sampledPoint, &sampledLight);

				float distance = (sampledPoint - ray.o).norm();
				ray.d = (sampledPoint - ray.o) / distance;
				yes = true;
				if (!scene.isOccluded(ray, distance)) {
					float solidAngle = abs(-ray.d.dot(sampledLight->n)) / nee_pdf / distance / distance; // rcp of the pdf
					if (solidAngle != 0.0f) {
						// TODO Can avoid 0.0f check by multiplying denom & numerator by solidAngle
						//std::cout << 1.0f / solidAngle / solidAngle / (1.0f / solidAngle / solidAngle + abs(ray.d.dot(n)) / M_PI * abs(ray.d.dot(n)) / M_PI) << '\n';
						color += sampledLight->material.emission * attenuation * material.albedo * abs(n.dot(ray.d)) / M_PI *
							1.0f / solidAngle / (1.0f / solidAngle / solidAngle + abs(ray.d.dot(n)) / M_PI * abs(ray.d.dot(n)) / M_PI);
						// balance heuristic with b = 2
						//exp(-density*distance); // MEDIUM TERM REMOVE LATER

					}
				}
			}

			//if Can Nee
			//	Sample point light
			//	Calc distance
			//	if !shadowcast & solidAngle not infinity
			//		calc solidangle
			//		color += emission * attenuation * albedo * brdf * lambert / brdf_pdf * mis_heuristic(1 / solidangle, brdf_pdf)
			//END NEE


			// TODO: Use onb instead
			Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0f, 0.0f, 1.0f));
			ray.d = rotMatrix * cosineSampleHemisphere(gen);
			float cos_t = std::max(ray.d.dot(n), 0.0f);
			if (yes) mis_brdf_pdf = cos_t / M_PI; else mis_brdf_pdf = -1.0f;
			attenuation = attenuation * material.albedo;/**cos_t/(cos_t/M_PI )*/ // * pi (Surface area) / (pi (lambertian albedo constant))
			break;
		}
		case specular: {
			// TODO Be able to handle materials with different refractive indexes
			float ior = 1.0f / 2.0f;
			float cos_t1 = -n.dot(ray.d);
			bool from_outside = cos_t1 > 0.0f;
			if (!from_outside) {
				// We're inside the specular object
				cos_t1 *= -1;
				ior = 1.0f / ior;
				n = -n;
			}

			// Check for total internal reflection, then choose to refract based on fresnel
			float cos_t2_sqr = 1.0f - ior * ior * (1.0f - cos_t1 * cos_t1);

			//TODO: See how much slower computing the real fresnel reflectance is. (will be twice as likely to call sqrt() than usual, but will be more accurate)
			if (cos_t2_sqr >= 0.0f &&
				gen.Uniform() > (/*(ior > 1.0f) ? schlickApprox(ior, cos_t1) :*/ schlickApprox(ior, from_outside ? cos_t1 : sqrt(cos_t2_sqr)))) {
				// Refraction through the specular surface
				float cos_t2 = sqrt(cos_t2_sqr);
				ray.d = ior * ray.d + (ior * cos_t1 - cos_t2) * n;
				ray.d /= ray.d.norm();
			}
			else {
				// Reflection off the specular surface
				ray.d = ray.d + (n * cos_t1 * 2);
			}
			attenuation = attenuation * material.albedo;
			mis_brdf_pdf = -1.0f;
			break;
		}
		}
		// Some diagnostic tools
		//i++;
		//c += last == obj && obj == light;
		//a += last == obj && obj == light;
		//last = obj;
	}
	//if (i > 1) cout << color << " " << i << endl;
	//if (a > 1) cout << color << " " << a << endl;
	//if (color.max() < 0.4) cout << color << endl;

	return color;
}
