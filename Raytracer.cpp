// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string.h>
#include <array>

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <random> // Replace with QMC generators

#include <chrono>
#include <ctime>  

#include <SFML/Graphics.hpp>

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"

#define WIDTH 1000
#define HEIGHT 1000
#define FOVX 120.0f
#define ACTIVE_SCENE 1

#define AMBIENT 0.3f
#define MAX_BOUNCES 7
#define SAMPLES_PER_PIXEL 80
#define BOUNCE_PROB 0.8f
#define HEMISPHERE_AREA (M_PI*2.0f)

/*
Things to try to remove noise
- Implement brightness based russian roulette
- better self intersection handling
- Investigate better light handling
  - Black body light? (basically just stop when the ray hits it)
  - Why does specular lighting cause noise?
  - Why does diffuse lighting make the light appear dark?
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

	Intersection castRay(Ray& ray) {
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

mt19937 gen;
uniform_real_distribution<float> pix(-0.5f, 0.5f);
uniform_real_distribution<float> uniform01(0.0f, 1.0f);
uniform_real_distribution<float> uniformAngle(0.0f, 2*M_PI);
Obj *light = NULL;
void buildScene(int i) {
	Material mirrorMat = { Vec3(1.0f), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(diffuse) };
	Material redMat = { Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) };
	Material greenMat = { Vec3(0.12f, 0.45f, 0.15f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };

	// Using specular light sauces creates a lot of noise
	Material lightMat = { Vec3(2.0f), 1.0f, Surface(reflective) };
	Material ovenMat = { Vec3(0.5f), 0.5f, Surface(diffuse) };
	switch (i) {
	default:
		break;
	case 0: // Plane
		scene.addObject(new Plane(Vec3(0.0f, 2.05f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
		//scene.addObject(new Sphere(Vec3(2.0f, 1.0f, -4.0f), 0.8f), lightMat);
		scene.addObject(new Sphere(Vec3(1.7f, 0.5f, -4.0f), 1.3f), specularMat);
		scene.addObject(new Box(Vec3(-1.4f, 0.5f, -3.0f), Vec3(-0.4f, 1.5f, -2.0f)), specularMat);
		break;
	case 1: // Cornell box
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), greenMat);
		scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), redMat);

		scene.addObject(new Box(Vec3(-600.0f, 400.0f, -600.0f), Vec3(-400.0f, 800.0f, -400.0f)), specularMat);
		scene.addObject(new Sphere(Vec3(400.0f, 440.0f, -600.0f), 200.3f), specularMat);
		//scene.addObject(new Box(Vec3(-200.0f, -800.0f, -500.0f), Vec3(200.0f, -750.0f, -100.0f)), diffuseMat);

		scene.addObject(new Box(Vec3(-200.0f, -800.0f, -500.0f), Vec3(200.0f, /*-795.0f*/ -799.9f, -100.0f)), lightMat);
		light = scene.objs.back();
		break;
	case 2: // Oven test
		/*
		* The oven test is any encosed room with surface emission 0.5 and reflectance 0.5. So we expect a pixel value
		* 0.5*(0.5 + 0.5(0.5 + 0.5(...)) = 1.
		*/
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), ovenMat);
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), ovenMat);
		scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), ovenMat);
		scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), ovenMat);
		scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), ovenMat);
		scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), ovenMat);
		break;
	}
	//scene.addObject(new Box(Vec3(1.0f, 1.0f, -3.0f), Vec3(1.0f, 1.0f, -3.0f) + Vec3(1.0f, 1.0f, 0.0f)), diffuseMat);
	/*

	/*for (Obj *obj : scene) {
		if (obj->emission > 0.0f) {
			lights.push_back(obj);
		}
	}*/
}

// Transform pixel coordinates to perspective rays
// Camera is at (0,0,0) facing (0,0,-1)
void cameraRay(float px, float py, Ray& ray) {
	float x = (2.0f * px - WIDTH) / WIDTH * tan(FOVX * M_PI / 180.0f / 2.0f);
	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	float y = (2.0f * py - HEIGHT) / HEIGHT * tan(((float) HEIGHT) / WIDTH * FOVX * M_PI / 180.0f / 2.0f);
	float z = -1.0f;
	ray.o = Vec3();
	ray.d = Vec3(x,y,z).normalized();

	//return Ray(Vec3(), Vec3(x, y, z));
}

Vec3 cosineSampleHemisphere() {
	float r_sqr = uniform01(gen);
	float phi = uniformAngle(gen);
	float r = sqrt(r_sqr);
	float x = cos(phi)*r;
	float y = sin(phi)*r;
	return Vec3(x, y, sqrt(1.0f - x*x - y*y));
}
Vec3 uniformSampleHemisphere() {
	float z = uniform01(gen);
	float theta = uniformAngle(gen);
	float r = sqrt(1.0f-z*z);
	return Vec3(cos(theta)*r, sin(theta)*r, z);
}

// cos_t is the dot product of the normal and the incident vectors
float schlickApprox(float r, float cos_t) {
	float R0 = (r - 1.0f) / (r + 1.0f);
	R0 *= R0;
	float x = 1.0f - cos_t;
	float x2 = x * x;
	return R0 + (1.0f - R0)*x2*x2*x;
}

int c = 0;
Vec3 rayTrace(Ray& ray, int depth) {
	Vec3 attenuation = Vec3(1.0f);
	Vec3 color = Vec3(0.0f);
	Obj *last = NULL;
	int i = 0;

	float mis_brdf_pdf = -1.0f; // Negative when mis was not used
	while (true) {
		Intersection intersection = scene.castRay(ray);
		
		//TODO: Weigh this by MIS
		if (isinf(intersection.t)) {
			float a = ray.d.dot(Vec3(0.2, -0.8, -0.4).normalized());
			if (a > 0.999) {

				color = color + attenuation*Vec3(5);
			} else if (a > 0.96) {
				color = color + attenuation*Vec3(5)*(a-0.96f)*(a - 0.96f)/ (0.999f - 0.96f) /(0.999f-0.96f);
			}
			color = color + attenuation*Vec3(0.5,0.70,0.8);
			break;
		}
		Vec3 hp = intersection.t*ray.d + ray.o;
		ray.o = hp;

		Obj *obj = intersection.hitObj;
		Vec3 n = obj->normal(hp);
		Vec3 emission = Vec3(obj->material.emission);

		if (emission.x > 0.0f) { // TODO USE vector max
			if (mis_brdf_pdf > 0.0f) {
				float solidAngle = abs(-ray.d.dot(obj->normal(hp))) * 400 * 400 / intersection.t / intersection.t;
				color = color + (obj->material.emission)*attenuation*
					mis_brdf_pdf * mis_brdf_pdf / (1 / solidAngle / solidAngle + mis_brdf_pdf * mis_brdf_pdf);
			} else {
				color = color + emission * attenuation;
			}
		}
		if (obj->material.surface == reflective) {
			float cos_t = ray.d.dot(n);

			ray.d = (ray.d - (n * cos_t * 2.0f)).normalized();
			attenuation = attenuation * obj->material.albedo;
			mis_brdf_pdf = -1.0f;
		} else if (obj->material.surface == diffuse) {
			float pl = 1.0f / 400.0f / 400.0f;
			if (light != obj) {
				Vec3 sampledLight = light->samplePoint(gen);
				ray.d = (sampledLight - ray.o).normalized();
				Intersection v = scene.castRay(ray);
				if ( light == v.hitObj) {
					float solidAngle = abs(-ray.d.dot(v.hitObj->normal(ray.o + v.t*ray.d)))/pl / v.t / v.t;

					//float pl = solidAngle * 1200.0f * 400.0f;
					// Gotta MIS this
					color = color + light->material.emission * attenuation * obj->material.albedo / M_PI * abs(n.dot(ray.d)) *
						1 / solidAngle /( 1/solidAngle/solidAngle + abs(ray.d.dot(n))/M_PI* abs(ray.d.dot(n)) / M_PI);
				}
			}
			//color = color + emission * attenuation;
			Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0f, 0.0f, 1.0f));
			ray.d = rotMatrix * cosineSampleHemisphere();
			float cos_t = abs(ray.d.dot(n));
			mis_brdf_pdf = cos_t / M_PI;

			//float solidAngle = 0.0;
			//intersection = scene.castRay(ray);
			//if (intersection.hitObj == light) {
			//	//return color;
			//	
			//	solidAngle = abs(-ray.d.dot(intersection.hitObj->normal(ray.o + intersection.t*ray.d))) * 400 * 400 / intersection.t / intersection.t;
			//	color = color + (intersection.hitObj->material.emission)*attenuation*obj->material.albedo*
			//		cos_t/M_PI * cos_t / M_PI /(1/solidAngle/solidAngle + cos_t/M_PI* cos_t / M_PI);
			//}
			attenuation = attenuation * obj->material.albedo;/**cos_t/(cos_t/M_PI )*/ // * pi (Surface area) / (pi (lambertian albedo constant))
		} else if (obj->material.surface == specular) {
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
				ray.d = (r*ray.d + (r*cos_t1 - cos_t2)*n).normalized();
			} else {
				// Reflection off the specular surface
				ray.d = (ray.d + (n * cos_t1 * 2)).normalized();
			}
			attenuation = attenuation * obj->material.albedo;
			mis_brdf_pdf = -1.0f;
		}

		// TODO: Move russian roulette
		//if (i > 1) {
		if (uniform01(gen) > BOUNCE_PROB) {
			break;
		}
		attenuation = attenuation / BOUNCE_PROB;
		i++;
		c += last == obj;
		last = obj;
	}
	return color;
}

void render(array<array<Vec3, WIDTH>, HEIGHT> *img) {
	auto t1 = chrono::high_resolution_clock::now();

	#pragma omp parallel for schedule(dynamic) num_threads(2)
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			Ray ray = Ray(Vec3(), Vec3());
			Vec3 colorVec = Vec3();
			for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
				cameraRay(px + pix(gen), py + pix(gen), ray);
				colorVec = colorVec + rayTrace(ray, MAX_BOUNCES);
			}
			(*img)[py][px] = colorVec / SAMPLES_PER_PIXEL;
		}
	}
	cout << c << endl;
	// Performance metrics
	float seconds = (chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0);
	cout << "Took " + to_string(seconds) + "s\n";
	cout << "Cast " << SAMPLES_PER_PIXEL * WIDTH*HEIGHT << " rays\n";
	cout << SAMPLES_PER_PIXEL * WIDTH*HEIGHT * (BOUNCE_PROB / (1.0f - BOUNCE_PROB) + 1) / seconds << " rays per second \n";
	cout << "Average bounces: " << BOUNCE_PROB/(1.0f - BOUNCE_PROB) << endl;
}

void toneMap(array<array<Vec3, WIDTH>, HEIGHT> *img) {
	//for (int py = 0; py < HEIGHT; py++) {
	//	for (int px = 0; px < WIDTH; px++) {
	//		(*img)[py][px] = Vec3(sqrt((*img)[py][px].x), sqrt((*img)[py][px].y), sqrt((*img)[py][px].z));
	//	}
	//}
}


string getDateTime() {
	// get the current time
	time_t t = time(0);

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

	// Create scene
	//cout << "Choose scene (0-" << MAX_SCENES << "): ";
	//int sceneIndex = 0;
	//cin >> sceneIndex;
	buildScene(ACTIVE_SCENE);

	array<array<Vec3, WIDTH>, HEIGHT>* img = new array<array<Vec3, WIDTH>, HEIGHT>();
	// Render scene
	render(img);



	// Tone map resulting image
	toneMap(img);

	sf::Image image;
	image.create(WIDTH, HEIGHT);
	// Load image
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			image.setPixel(px, py, (*img)[py][px].tosRGB());
		}
	}

	// Init texture
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer");
	sf::Texture texture;
	sf::Sprite sprite;
	texture.loadFromImage(image);
	sprite.setTexture(texture);

	ofstream ofile;
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
					cout << "Saving render" << endl;

					image.saveToFile(".\\renders\\render" + getDateTime() + ".png");
					ofile.open(".\\renders\\renderData" + getDateTime() + ".txt");
					for (int py = 0; py < HEIGHT; py++) {
						for (int px = 0; px < WIDTH; px++) {
							Vec3 v = (*img)[py][px];
							ofile << (int)v.x << ',' << (int)v.y << ',' << (int)v.z << ' ';
						}
						ofile << endl;
					}
					ofile.close();

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
	delete img;
	return 0;
}