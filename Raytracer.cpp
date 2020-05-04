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


#define WIDTH 1000
#define HEIGHT 1000
#define FOVX (50.0/180.0*M_PI)

using namespace std;

vector<Obj *> scene;
vector<Obj *> lights;

void buildScene() {
	scene.push_back(new Sphere(Vec3(0.0, 0.0, -4.0), 2.0));

	for (Obj *obj : scene) {
		if (obj->emission > 0.0) {
			lights.push_back(obj);
		}
	}
}

// Transform pixel coordinates to perspective rays
Ray cameraRay(int px, int py) {
	double x = (2.0 * px - WIDTH) / WIDTH * tan(FOVX);
	double y = (2.0 * py - HEIGHT) / HEIGHT * tan(((double) HEIGHT) / WIDTH * FOVX);
	double z = -1;
	return Ray(Vec3(), Vec3(x, y, z));
}

// TODO Get rid of tuples
tuple<double, Obj *> castRay(Ray ray) {
	double minT = -1.0;
	Obj *hitObj = NULL;
	for (Obj *obj : scene) {
		double t = obj->intersect(ray);
		if (t > 0.0 && (minT < 0.0 || t < minT)) {
			minT = t;
			hitObj = obj;
		}
	}
	return make_tuple(minT, hitObj);
}

void draw(sf::Image &image) {
	auto t1 = chrono::high_resolution_clock::now();

	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			tuple<double, Obj *> rayHit = castRay(cameraRay(px, py));
			if (get<0>(rayHit) > 0.0) {
				image.setPixel(px, py, sf::Color(255, 255, 255));
			} else {
				image.setPixel(px, py, sf::Color(0, 0, 0));
			}
		}
	}

	cout << "Took " + to_string(chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0) + "s\n";
}

int main() {
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer");
	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;
	image.create(WIDTH, HEIGHT);
	buildScene();
	draw(image);
	texture.loadFromImage(image);
	sprite.setTexture(texture);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		window.clear();

		window.draw(sprite);


		window.display();
	}
	
	return 0;
}
