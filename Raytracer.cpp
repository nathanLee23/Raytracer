// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <random> // Replace with QMC generators

#include <chrono>
#include <ctime>  

#include <SFML/Graphics.hpp>

#include "Vec3.h"
#include "Material.h"
#include "Obj.h"
#include "Sphere.h"
#include "Plane.h"
#include "constants.h"
#include "Matrix3.h"

#define WIDTH 1000
#define HEIGHT 1000
#define FOVX (50.0f/180.0f*M_PI)

#define AMBIENT 300.0f
#define MIN_BOUNCES 7
#define SAMPLES_PER_PIXEL 50

/*
Easy
Implement AABB
Try Vec4
Russian roulette
Try using plane reflection to generate rotated hemisphere over rotation matrices. Branching can be avoided by being clever
Try precomputing the rotation matrix since we have that the second Vec is always (0,0f,1)

Medium
The BVH is a static tree which might be optimizable.
The BVH being static could be made into an array which is likely faster. See PBRT book
- Shadow rays, operate on a line segment, if ANY intersection is found between the 2 points then that is sufficient and you can early terminate,
- Light rays must find the closest point of interesction so all points must be investigated.
Research self intersection in raytracing gems

Hard
See jacco blog for kernel optimization (use small kernels over a singular megakernel).
Data oriented design (DOD), means that we should not use OOP and instead store each type of 'object' in their own array which reduces branching and improves
memory locality. This means all materials and geometric primitives should be processed in their own array.
Could use CUDA SIMD intrinsics for even more speed.
Metropolis-Hastings
Multiple importance sampling
*/

using namespace std;

struct Intersection {
	float t;
	Obj * hitObj;
	Material material;
};

class Scene {
public:
	vector<Obj *> objs;

	vector<Obj *> lights;

	void addObject(Obj *obj, Material material) {
		obj->material = material;
		objs.push_back(obj);
	}

	Intersection castRay(Ray ray) {
		float minT = INFINITY;
		Obj *hitObj = NULL;
		for (Obj *obj : objs) {
			float t = obj->intersect(ray);
			if (t > EPS && t < minT) {
				minT = t;
				hitObj = obj;
			}
		}
		return {minT, hitObj};
	}

};

Scene scene;

default_random_engine generator;
uniform_real_distribution<float> pix(-0.5f, 0.5f); // Need a better name for this
uniform_real_distribution<float> hemisphere(0.0f, 1.0f); // Need a better name for this

void buildScene() {
	Material mirrorMat = { Vec3(), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.9f,0.9f,0.9f), 0.0f, Surface(diffuse) };
	Material redMat = { Vec3(M_PI / 2.0f,0.3f,0.3f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };
	Material lightMat = { Vec3(), 320.0f, Surface(reflective) };

	scene.addObject(new Plane(Vec3(0.0f, 2.05f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
	scene.addObject(new Sphere(Vec3(0.0f, 0.5f, -4.0f), 1.3f), specularMat);
	scene.addObject(new Sphere(Vec3(1.7f, 0.5f, -6.0f), 1.3f), diffuseMat);
	/*
	scene.addObject(new Plane(Vec3(0.0f, 0.0f, -7.0f), Vec3(0.0f, 0.0f, 1.0)), diffuseMat);
	scene.addObject(new Plane(Vec3(0.0f, 0.0f, 2.0f), Vec3(0.0f, 0.0f, -1.0)), diffuseMat);*//*
	scene.addObject(new Plane(Vec3(0.0f, -4.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), lightMat);
	scene.addObject(new Plane(Vec3(4.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), diffuseMat);
	scene.addObject(new Plane(Vec3(-4.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), redMat);*/

	/*for (Obj *obj : scene) {
		if (obj->emission > 0.0f) {
			lights.push_back(obj);
		}
	}*/
}

// Transform pixel coordinates to perspective rays
// Camera is at (0,0f,0) facing (0,0f,-1)
Ray cameraRay(float px, float py) {
	float x = (2.0f * px - WIDTH) / WIDTH * tan(FOVX);
	float y = (2.0f * py - HEIGHT) / HEIGHT * tan(((float) HEIGHT) / WIDTH * FOVX);
	float z = -1.0f;
	return Ray(Vec3(), Vec3(x, y, z));
}

Vec3 sampleHemisphere() {
	float u1 = hemisphere(generator);
	float u2 = hemisphere(generator);
	float r = sqrt(1.0f - u1 * u1);
	float phi = u2 * M_PI * 2.0f;
	return Vec3(cos(phi)*r, sin(phi)*r, u1);
}

float shlickApprox(float r, float R0, float cos_t2) {
	float x = 1.0f - cos_t2;
	float x2 = x * x;
	return R0 + (1.0f - R0)*x2*x2*x;
}

Vec3 rayTrace(Ray ray, int depth) {
	if (depth <= 0) {
		return Vec3();
	}
	Intersection intersection = scene.castRay(ray);
	if (isinf(intersection.t)) {
		return Vec3(AMBIENT);
	}

	Vec3 hp = intersection.t*ray.d + ray.o;
	ray.o = hp;

	Obj *obj = intersection.hitObj;
	Vec3 n = obj->normal(hp);


	Vec3 clr = Vec3(obj->material.emission);

	if (obj->material.surface == reflective) {
		float cos_t = ray.d.dot(n);

		ray.d = (ray.d - (n * cos_t * 2.0f)).normalized();
		return clr + rayTrace(ray, depth - 1);
	} else if (obj->material.surface == diffuse) {
		// This line is optimizable since the vector to rotate is constant
		Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0f, 0.0f, 1.0f));
		ray.d = rotMatrix*sampleHemisphere();
		float cos_t = ray.d.dot(n);
		return clr + rayTrace(ray, depth - 1)*obj->material.albedo * cos_t;
	} else if (obj->material.surface == specular) {
		float r = 1.0f / 1.3f;
		float cos_t1 = -n.dot(ray.d);
		if (cos_t1 < 0.0f) {
			// We're inside the specular object
			r = 1 / r;
			cos_t1 *= -1;
			n = -n;
		}
		float R0 = (1.0f - 1.3f) / (1.0f + 1.3f);
		R0 *= R0;
		float cos_t2 = sqrt(1.0f - r * r*(1.0f - cos_t1 * cos_t1));
		float R = shlickApprox(r, R0, cos_t2);
		if (cos_t2 >= 0.0f && hemisphere(generator) > R) {
			ray.d = (r*ray.d + (r*cos_t1 - cos_t2)*n).normalized();
		} else {
			ray.d = (ray.d - (n * cos_t1 * 2)).normalized();
		}
		return clr + rayTrace(ray, depth - 1) * obj->material.albedo;
	}
	return clr;
}

void draw(sf::Image &image) {
	auto t1 = chrono::high_resolution_clock::now();

	#pragma omp parallel for schedule(dynamic) num_threads(2)
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			Vec3 colorVec = Vec3();
			for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
				colorVec = colorVec + rayTrace(cameraRay(px + pix(generator), py + pix(generator)), MIN_BOUNCES);
			}
			image.setPixel(px, py, (colorVec/SAMPLES_PER_PIXEL).toColor()); // Not threadsafe
		}
	}

	cout << "Took " + to_string(chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0) + "s\n";
}

string getDateTime() {
	time_t t = time(0);   // get time now
	tm now;
	localtime_s(&now, &t);
	string s = to_string(now.tm_year + 1900)
		+ to_string(now.tm_mon + 1)
		+ to_string(now.tm_mday)
		+ to_string(now.tm_hour)
		+ to_string(now.tm_min)
		+ to_string(now.tm_sec);
	return s;
}

int main() {
	if (!std::numeric_limits<float>::is_iec559) {
		cout << "Machine architecture must implement IEEE 754.\n";
		return 0;
	}
	// Init texture
	sf::Image image;
	image.create(WIDTH, HEIGHT);

	// Render scene
	buildScene();
	draw(image);

	// Load image
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer");
	sf::Texture texture;
	sf::Sprite sprite;
	texture.loadFromImage(image);
	sprite.setTexture(texture);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				switch (event.key.code) {
				case sf::Keyboard::S:
					if (image.saveToFile(".\\renders\\render" + getDateTime() + ".png"))
						cout << "Render saved" << endl;
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
	}
	
	return 0;
}