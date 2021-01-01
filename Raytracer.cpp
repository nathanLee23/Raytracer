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
#include <random> // Replace with QMC generators

#include <chrono>
#include <ctime>  

#include <omp.h>
#include <SFML/Graphics.hpp>

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"

#define WIDTH 800
#define HEIGHT 800
#define FOVX 120.0f
#define ACTIVE_SCENE 3

#define AMBIENT 0.3f
#define SAMPLES_PER_PIXEL 1
#define THREADS 2

// Unused
#define MIN_BOUNCES 2
#define BOUNCE_PROB 0.8f
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
		(*pixels)[py][px] = ((*pixels)[py][px] * (float) (sample_count) + color) / ((float) (sample_count + 1));
	}
};
struct Intersection {
	float t;
	Obj * hitObj;
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
		Obj *hitObj = nullptr;
		for (Obj *obj : objs) {
			float t = obj->intersect(ray);
			if (t > EPS && t < minT) {
				minT = t;
				hitObj = obj;
			}
		}
		return {minT, hitObj};
	}
/*
	float mesh_intersect(Ray &ray, Vec3 &n, Obj **hitObj) {

	}*/
};

Scene scene;

uniform_real_distribution<float> pix(-0.5f, 0.5f);
uniform_real_distribution<float> uniform01(0.0f, 1.0f);
uniform_real_distribution<float> uniformAngle(0.0f, 2.0f*M_PI);
Obj *light = nullptr;
void buildScene(int i) {
	Material mirrorMat = { Vec3(1.0f), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(diffuse) };
	Material redMat = { Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) };
	Material greenMat = { Vec3(0.12f, 0.45f, 0.15f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };

	// Using specular light sauces creates a lot of noise
	Material lightMat = { Vec3(1.0f), 16.0f, Surface(specular) };
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

		scene.addObject(new Sphere(Vec3(400.0f, 440.0f, -600.0f), 200.3f), specularMat);
		scene.addObject(new Box(Vec3(-200.0f, -800.0f, -500.0f), Vec3(200.0f, -750.0f, -100.0f)), specularMat);

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
	case 3:
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, -800.0f), Vec3(0.0f, 0.0f, 1.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, 0.0f, 800.0f), Vec3(0.0f, 0.0f, -1.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, 800.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f)), diffuseMat);
		scene.addObject(new Plane(Vec3(0.0f, -800.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0)), diffuseMat);
		scene.addObject(new Plane(Vec3(800.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0)), greenMat);
		scene.addObject(new Plane(Vec3(-800.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0)), redMat);

		scene.addObject(new Sphere(Vec3(0.0f, -040.0f, -500.0f), 200.3f), diffuseMat);
		scene.addObject(new Box(Vec3(-50.0f, -800.0f, -650.0f), Vec3(50.0f, /*-795.0f*/ -799.9f, -550.0f)), lightMat);
		light = scene.objs.back();
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
	float x = (2.0f * px - WIDTH) / WIDTH * tanf(FOVX * M_PI / 180.0f / 2.0f);
	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	float y = (2.0f * py - HEIGHT) / HEIGHT * tanf(((float) HEIGHT) / WIDTH * FOVX * M_PI / 180.0f / 2.0f);
	float z = -1.0f;
	ray.o = Vec3(0.0,0.0f,0.0f);
	ray.d = Vec3(x,y,z).normalized();

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
	return 1/4.0f/M_PI*(1.0f-g*g)/(1.0f+g*g+2.0f*g*powf(cos_t,1.5f));
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
	int a = 0;
	float mis_brdf_pdf = -1.0f; // Negative when mis was not used
	while (true) {
		Intersection intersection = scene.castRay(ray);
		
		//TODO: Weigh this by MIS
		if (isinf(intersection.t)) {
			float a = ray.d.dot(Vec3(0.2f, -0.8f, -0.4f).normalized());
			if (a > 0.999) {

				color = color + attenuation*Vec3(5);
			} else if (a > 0.96) {
				color = color + attenuation*Vec3(5)*(a-0.96f)*(a - 0.96f)/ (0.999f - 0.96f) /(0.999f-0.96f);
			}
			color = color + attenuation*Vec3(0.5f,0.70f,0.8f);
			break;
		}

		//float density = 0.0005f;
		float density = 0.0f;

		Vec3 hp = intersection.t*ray.d + ray.o;

		Obj *obj = intersection.hitObj;
		Vec3 n = obj->normal(hp);
		Vec3 emission = Vec3(obj->material.emission);
		if (emission.max() > 0.0f) {
			if (mis_brdf_pdf > 0.0f) {
				float solidAngle = abs(-ray.d.dot(obj->normal(hp))) * 100 * 100 / intersection.t / intersection.t;
				color = color + (obj->material.emission)*attenuation *
					mis_brdf_pdf * mis_brdf_pdf / (1 / solidAngle / solidAngle + mis_brdf_pdf * mis_brdf_pdf) *
					exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER
				
			} else {
				color = color + emission * attenuation
					* exp(-density * intersection.t); // MEDIUM TERM REMOVE LATER
			}
		}
		//float medium_t = -log(1.0f - uniform01(gen)) / density;
		float medium_t = INFINITY;
		if (intersection.t > medium_t) {
			//hit medium

			ray.o = medium_t * ray.d + ray.o;

			// MIS direct light
			Vec3 sampledLight = light->samplePoint(gen);
			ray.d = (sampledLight - ray.o).normalized();
			Intersection v = scene.castRay(ray);
			if (light == v.hitObj) {
				float pl = 1.0f / 100.0f / 100.0f;
				float solidAngle = abs(-ray.d.dot(v.hitObj->normal(ray.o + v.t*ray.d))) / pl / v.t / v.t;

				color = color + light->material.emission * attenuation * obj->material.albedo / M_PI * abs(n.dot(ray.d)) *
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
		float bounce_probability = min(attenuation.max(), 0.95f);
		if (uniform01(gen) > bounce_probability) {
			break;
		}
		attenuation = attenuation / bounce_probability;

		if (obj->material.surface == reflective) {
			float cos_t = ray.d.dot(n);

			ray.d = (ray.d - (n * cos_t * 2.0f));
			attenuation = attenuation * obj->material.albedo;
			mis_brdf_pdf = -1.0f;
		} else if (obj->material.surface == diffuse) {
			// TODO actually do this properly
			float pl = 1.0f / 100.0f / 100.0f;
			if (light != obj) {
				Vec3 sampledLight = light->samplePoint(gen);
				ray.d = (sampledLight - ray.o).normalized();
				Intersection v = scene.castRay(ray);
				if ( light == v.hitObj) {
					float solidAngle = abs(-ray.d.dot(v.hitObj->normal(ray.o + v.t*ray.d)))/pl / v.t / v.t;

					color = color + light->material.emission * attenuation * obj->material.albedo / M_PI * abs(n.dot(ray.d)) *
						1 / solidAngle /( 1/solidAngle/solidAngle + abs(ray.d.dot(n))/M_PI* abs(ray.d.dot(n)) / M_PI) *
						exp(-density*v.t); // MEDIUM TERM REMOVE LATER
				}
			}
			//color = color + emission * attenuation;
			Matrix3 rotMatrix = rotMatrixVectors(n, Vec3(0.0f, 0.0f, 1.0f));
			ray.d = rotMatrix * cosineSampleHemisphere(gen);
			float cos_t = abs(ray.d.dot(n));
			mis_brdf_pdf = cos_t / M_PI;
			attenuation = attenuation * obj->material.albedo;/**cos_t/(cos_t/M_PI )*/ // * pi (Surface area) / (pi (lambertian albedo constant))
		} else if (obj->material.surface == specular) {
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
			attenuation = attenuation * obj->material.albedo;
			mis_brdf_pdf = -1.0f;
		}

		// Some diagnostic tools
		//i++;
		//c += last == obj && obj == light;
		//a += last == obj && obj == light;
		//last = obj;
	}
	//if (a > 1) cout << color << " " << a << endl;
	return color;
}

random_device rd;
mt19937 gens[THREADS];
auto t1 = chrono::high_resolution_clock::now();
void render(Img &img) {
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
	cout << "Took " + to_string(seconds) + "s\n";
	/*cout << "Cast " << SAMPLES_PER_PIXEL * WIDTH*HEIGHT << " rays\n";
	cout << SAMPLES_PER_PIXEL * WIDTH*HEIGHT / seconds << " samples per second \n";*/
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
	string s = to_string(now.tm_year + 1900)
		+ to_string(now.tm_mon + 1)
		+ to_string(now.tm_mday)
		+ to_string(now.tm_hour)
		+ to_string(now.tm_min)
		+ to_string(now.tm_sec);
	return s;
}

void render_loop(Img &img) {
	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < THREADS; i++) gens[i] = mt19937(hash<int>{}(i));
	while (true) {
		render(img);
		img.sample_count++;
		cout << "Accumulated " << SAMPLES_PER_PIXEL*img.sample_count << " samples per pixel" << endl;
	}
}

void process_image(sf::Image &image, array<array<Vec3, WIDTH>, HEIGHT> *img) {
	toneMap(img);
	// Load image
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			image.setPixel(px, py, (*img)[py][px].tosRGB());
		}
	}
}

void gui_thread(Img &img) {
	// Tone map resulting image

	sf::Image image;
	image.create(WIDTH, HEIGHT);

	//process_image(image, img);

	// Init texture
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer");
	sf::Texture texture;
	sf::Sprite sprite;


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
					cout << "Saving render" << endl;

					image.saveToFile(".\\renders\\render" + getDateTime() + ".png");
					ofile.open(".\\renders\\renderData" + getDateTime() + ".txt");
					for (int py = 0; py < HEIGHT; py++) {
						for (int px = 0; px < WIDTH; px++) {
							Vec3 v = (*img.pixels)[py][px];
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

		process_image(image, img.pixels);
		texture.loadFromImage(image);
		sprite.setTexture(texture);
	}
}

int main() {
	if (!std::numeric_limits<float>::is_iec559) {
		cout << "Machine architecture must implement IEEE 754.\n";
		return 0;
	}
	auto t1 = chrono::high_resolution_clock::now();
	load_mesh("geometry/teapot.obj");
	float seconds = (chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0);
	cout << "Took " + to_string(seconds) + "s\n";
	// Create scene
	//cout << "Choose scene (0-" << MAX_SCENES << "): ";
	//int sceneIndex = 0;
	//cin >> sceneIndex;
	buildScene(ACTIVE_SCENE);

	Img img = Img();
	img.clear();

	// Render scene
	thread render_thread = thread(render_loop, ref(img));
	gui_thread(img);
	exit(0);
	return 0;
}