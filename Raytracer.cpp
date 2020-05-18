// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>  
#include <string.h>
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vec3.h"
#include "Material.h"
#include "Obj.h"
#include "Sphere.h"
#include "Plane.h"
#include "constants.h"
#include "Matrix3.h"

#include <random> // Replace with QMC generators

#define WIDTH 1000
#define HEIGHT 1000
#define FOVX (50.0/180.0*M_PI)

#define AMBIENT 300.0
#define MIN_BOUNCES 7
#define SAMPLES 40

using namespace std;

struct Intersection {
	double t;
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
		double minT = -1.0;
		Obj *hitObj = NULL;
		for (Obj *obj : objs) {
			double t = obj->intersect(ray);
			if (t > EPS && (minT < 0.0 || t < minT)) {
				minT = t;
				hitObj = obj;
			}
		}
		return {minT, hitObj};
	}

};

Scene scene;

default_random_engine generator;
uniform_real_distribution<double> pix(-0.5, 0.5); // Need a better name for this
uniform_real_distribution<double> hemisphere(0.0, 1.0); // Need a better name for this

void buildScene() {
	Material mirrorMat = { Vec3(), 0.0, Surface(reflective) };
	Material diffuseMat = { Vec3(9.0,9.0,9.0), 0.0, Surface(diffuse) };
	Material specularMat = { Vec3(), 0.0, Surface(specular) };
	Material lightMat = { Vec3(), 320.0, Surface(reflective) };
	scene.addObject(new Sphere(Vec3(0, 0.5, -4.0), 0.3), specularMat);
	scene.addObject(new Sphere(Vec3(0, 0.5, -4.0), 0.6), specularMat);
	scene.addObject(new Sphere(Vec3(0, 0.5, -4.0), 1.3), specularMat);
	//scene.addObject(new Plane(Vec3(0.0, 0.0, -7.0), Vec3(0.0, 0.0, 1.0)), diffuseMat);
	//scene.addObject(new Plane(Vec3(0.0, 0.0, 2.0), Vec3(0.0, 0.0, -1.0)), diffuseMat);
	scene.addObject(new Plane(Vec3(0.0, 2.05, 0.0), Vec3(0.0, -1.0, 0.0)), diffuseMat);
	//scene.addObject(new Plane(Vec3(0.0, -4.0, 0.0), Vec3(0.0, 1.0, 0.0)), lightMat);
	//scene.addObject(new Plane(Vec3(4.0, 0.0, 0.0), Vec3(-1.0, 0.0, 0.0)), diffuseMat);
	//scene.addObject(new Plane(Vec3(-4.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0)), diffuseMat);

	/*for (Obj *obj : scene) {
		if (obj->emission > 0.0) {
			lights.push_back(obj);
		}
	}*/
}

// Transform pixel coordinates to perspective rays
// Camera is at (0,0,0) facing (0,0,-1)
Ray cameraRay(double px, double py) {
	double x = (2.0 * px - WIDTH) / WIDTH * tan(FOVX);
	double y = (2.0 * py - HEIGHT) / HEIGHT * tan(((double) HEIGHT) / WIDTH * FOVX);
	double z = -1;
	return Ray(Vec3(), Vec3(x, y, z));
}

Vec3 sampleHemisphere() {
	double u1 = hemisphere(generator);
	double u2 = hemisphere(generator);
	double r = sqrt(1.0 - u1 * u1);
	double phi = u2 * M_PI * 2.0;
	return Vec3(cos(phi)*r, sin(phi)*r, u1);
}

double shlickApprox(double r, double R0, double cos_t2) {
	double x = 1 - cos_t2;
	double x2 = x * x;
	return R0 + (1 - R0)*x2*x2*x;
}

Vec3 rayTrace(Ray ray, int depth) {
	if (depth <= 0) {
		return Vec3();
	}
	Intersection intersection = scene.castRay(ray);
	if (intersection.t < 0.0 ) {
		return Vec3(AMBIENT);
	}
	Vec3 hp = intersection.t*ray.d + ray.o;
	ray.o = hp;

	Obj *obj = intersection.hitObj;
	Vec3 n = obj->normal(hp);


	Vec3 clr = Vec3(obj->material.emission);

	if (obj->material.surface == reflective) {
		double cos_t = ray.d.dot(n);

		ray.d = (ray.d - (n * cos_t * 2)).normalized();
		return clr + rayTrace(ray, depth - 1);
	} else if (obj->material.surface == diffuse) {
		Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0, 0.0, 1.0));
		ray.d = rotMatrix*sampleHemisphere();
		double cos_t = ray.d.dot(n);
		return clr + rayTrace(ray, depth - 1)*obj->material.albedo * 0.1 * cos_t;
	} else if (obj->material.surface == specular) {
		double r = 1.0 / 1.3;
		double cos_t1 = -n.dot(ray.d);
		if (cos_t1 < 0.0) {
			// We're inside the specular object
			r = 1 / r;
			cos_t1 *= -1;
			n = -n;
		}
		double R0 = (1.0 - 1.3) / (1.0 + 1.3);
		R0 *= R0;
		double cos_t2 = sqrt(1 - r * r*(1 - cos_t1 * cos_t1));
		double R = shlickApprox(r, R0, cos_t2);
		if (cos_t2 >= 0 && hemisphere(generator) > R) {
			ray.d = (r*ray.d + (r*cos_t1 - cos_t2)*n).normalized();
		} else {
			ray.d = (ray.d - (n * cos_t1 * 2)).normalized();
		}
		return clr + rayTrace(ray, depth - 1);
	}
	return clr;
}

void draw(sf::Image &image) {
	auto t1 = chrono::high_resolution_clock::now();

	#pragma omp parallel for schedule(dynamic) num_threads(2)
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			Vec3 colorVec = Vec3();
			for (int sample = 0; sample < SAMPLES; sample++) {
				colorVec = colorVec + rayTrace(cameraRay(px + pix(generator), py + pix(generator)), MIN_BOUNCES);
			}
			image.setPixel(px, py, (colorVec/SAMPLES).toColor()); // Not threadsafe
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