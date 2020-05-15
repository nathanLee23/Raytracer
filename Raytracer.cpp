// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <string.h>
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vec3.h"
#include "Obj.h"
#include "Sphere.h"
#include "Plane.h"
#include "constants.h"

#include <random> // Replace with QMC generators

#define WIDTH 1000
#define HEIGHT 1000
#define FOVX (50.0/180.0*M_PI)

#define MAX_BOUNCES 8
#define SAMPLES 5

using namespace std;

struct Intersection {
	double t;
	Obj * hitObj;
};

vector<Obj *> scene;
vector<Obj *> lights;

default_random_engine generator;
uniform_real_distribution<double> pix(-0.5, 0.5); // Need a better name for this

void buildScene() {
	scene.push_back(new Sphere(Vec3(0.0, 0.0, -4.0), 1.5));
	scene.push_back(new Plane(Vec3(0.0, 0.0, -7.0), Vec3(0.0, 0.0, 1.0)));
	scene.push_back(new Plane(Vec3(0.0, 0.0, 2.0), Vec3(0.0, 0.0, -1.0)));
	scene.push_back(new Plane(Vec3(0.0, 4.0, 0.0), Vec3(0.0, -1.0, 0.0)));
	scene.push_back(new Plane(Vec3(0.0, -4.0, 0.0), Vec3(0.0, 1.0, 0.0)));
	scene.push_back(new Plane(Vec3(4.0, 0.0, 0.0), Vec3(-1.0, 0.0, 0.0)));
	scene.push_back(new Plane(Vec3(-4.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0)));

	/*for (Obj *obj : scene) {
		if (obj->emission > 0.0) {
			lights.push_back(obj);
		}
	}*/
}

// Transform pixel coordinates to perspective rays
// Camera is at (0,0) facing (0,0,-1)
Ray cameraRay(double px, double py) {
	double x = (2.0 * px - WIDTH) / WIDTH * tan(FOVX);
	double y = (2.0 * py - HEIGHT) / HEIGHT * tan(((double) HEIGHT) / WIDTH * FOVX);
	double z = -1;
	return Ray(Vec3(), Vec3(x, y, z));
}

Intersection castRay(Ray ray) {
	double minT = -1.0;
	Obj *hitObj = NULL;
	for (Obj *obj : scene) {
		double t = obj->intersect(ray);
		if (t > EPS && (minT < 0.0 || t < minT)) {
			minT = t;
			hitObj = obj;
		}
	}
	return { minT, hitObj }; // Change to alter arguments
}


Vec3 rayTrace(Ray ray, int depth) {
	Vec3 clr = Vec3();
	if (depth <= 0) {
		return clr;
	}
	Intersection intersection = castRay(ray);
	if (intersection.t < EPS) {
		return clr;
	}
	Vec3 hp = intersection.t*ray.d + ray.o;
	Vec3 n = intersection.hitObj->normal(hp);
	clr = Vec3(intersection.hitObj->emission);

	double cost = ray.d.dot(n);

	ray.o = hp;
	ray.d = (ray.d - (n * cost * 2)).normalized();
	return clr + rayTrace(ray, depth - 1);
}

void draw(sf::Image &image) {
	auto t1 = chrono::high_resolution_clock::now();

	//#pragma omp parallel for
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			Vec3 colorVec = Vec3();
			for (int sample = 0; sample < SAMPLES; sample++) {
				colorVec = colorVec + rayTrace(cameraRay(px + pix(generator), py + pix(generator)), MAX_BOUNCES);
			}
			image.setPixel(px, py, (colorVec/SAMPLES).toColor()); // Not threadsafe
		}
	}

	cout << "Took " + to_string(chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0) + "s\n";
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
					if (image.saveToFile("render.png"))
						cout << "Render saved" << endl;
					else
						cout << "Failed to save render" << endl;
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
