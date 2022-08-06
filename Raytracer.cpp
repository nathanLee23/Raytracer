// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <thread>
#include <functional>
#include <string.h>
#include <atomic>
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
#include "Camera.h"
#include "Scene.h"
#include "Integrator.h"

#define FOVX 60.0f
#define ACTIVE_SCENE 3

#define AMBIENT 0.3f
#define SAMPLES_PER_PIXEL 4
#define THREADS 2

using namespace std::chrono_literals;

std::atomic<bool> should_reset;

class Img {
public:
	std::array<std::array<Vec3, WIDTH>, HEIGHT> &pixels; // TODO replace with something else, reference member is bad practice
	unsigned int sample_count;
	Img() : pixels(*(new std::array<std::array<Vec3, WIDTH>, HEIGHT>())) {
		clear();
	}
	
	void clear() {
		for (int py = 0; py < HEIGHT; py++) {
			for (int px = 0; px < WIDTH; px++) {
				pixels[py][px] = Vec3();
			}
		}
		sample_count = 0;
	}

	void update(int px, int py, Vec3 color) {
		// TODO use Kahan summation?
		pixels[py][px] = (pixels[py][px] * static_cast<float>(sample_count) + color) / static_cast<float>(sample_count + 1);
	}
};
auto camera = PerspectiveCamera(FOVX);
Scene scene = Scene();

void buildScene(int i) {
	Material mirrorMat = { Vec3(1.0f), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(diffuse) };
	Material varnishMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(varnish) };
	Material redMat = { Vec3(0.65f, 0.05f, 0.05f), 0.0f, Surface(diffuse) };
	Material greenMat = { Vec3(0.12f, 0.45f, 0.15f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };

	// Using specular light sauces creates a lot of noise
	Material lightMat = { Vec3(1.0f), 16.0f, Surface(specular) };
	Material ovenMat = { Vec3(0.5f), 0.5f, Surface(diffuse) };

	scene.addMesh(scene.load_mesh("geometry/CornellBox-Original.obj"));
	scene.meshes[0]->materials[7].emission = 4.0f;
	scene.meshes[0]->materials[6] = mirrorMat;

#ifdef SPHERES
	scene.spheres.push_back(Sphere(Vec3(0.0f,1.3f,0.5f), 0.3f));
	scene.spheres[0].material = specularMat;
#endif

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

	scene.commit();
}

std::random_device rd;
Pcg gens[THREADS];

template<typename I>
void render(Img &img, I&& integrator) {
	auto t1 = std::chrono::high_resolution_clock::now();
	if (img.sample_count == 0) {
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int py = 0; py < HEIGHT; py++) {
			Pcg& gen = gens[omp_get_thread_num()];
			for (int px = 0; px < WIDTH; px++) {
				Ray ray = Ray(Vec3(), Vec3());
				Vec3 colorVec = Vec3();
				for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
					camera.pixelToCameraRay(px + gen.Uniform()-0.5f, py + gen.Uniform()-0.5f, ray, gen);
					colorVec += integrator(scene, ray, gen);
				}
				img.pixels[py][px] = colorVec/SAMPLES_PER_PIXEL;
			}
		}
	} else {
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int py = 0; py < HEIGHT; py++) {
#ifdef RT_DEBUG
			// Enable exception on 'bad' float operations
			unsigned control_word;
			_controlfp_s(&control_word, ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_UNDERFLOW),
				_MCW_EM);
#endif
			Pcg& gen = gens[omp_get_thread_num()];
			for (int px = 0; px < WIDTH; px++) {
				Ray ray = Ray(Vec3(), Vec3());
				Vec3 colorVec = Vec3();
				for (int sample = 0; sample < SAMPLES_PER_PIXEL; sample++) {
					camera.pixelToCameraRay(px + gen.Uniform()-0.5f, py + gen.Uniform()-0.5f, ray, gen);
					colorVec += integrator(scene, ray, gen);
				}
				img.update(px, py, colorVec/SAMPLES_PER_PIXEL);
			}
		}

	}
		
	// Performance metrics
	float seconds = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count() / 1000000.0;
	printf("Took %f seconds\n", seconds);
	printf("Cast %d rays\n", SAMPLES_PER_PIXEL * WIDTH*HEIGHT);
	printf("%f samples per second \n", SAMPLES_PER_PIXEL * WIDTH * HEIGHT / seconds);
	printf("%f ms per sample \n", seconds * 1'000'000.0f / (SAMPLES_PER_PIXEL * WIDTH * HEIGHT));
}

void toneMap(std::array<std::array<Vec3, WIDTH>, HEIGHT> &img) {
	//for (int py = 0; py < HEIGHT; py++) {
	//	for (int px = 0; px < WIDTH; px++) {
	//		img[py][px] = Vec3(img[py][px].x, img[py][px].y, img[py][px].z);
	//	}
	//}
}

std::string getDateTime() {
	// get the current time
	time_t t = time(0);

	tm now;
	localtime_s(&now, &t);
	char buffer[15];
	std::strftime(buffer, 32, "%Y%m%d%H%M%S", &now);
	return std::string(buffer);
}

void render_loop(Img &img) {
	while (true) {
		for (int i = 0; i < THREADS; i++) {
			gens[i] = Pcg(std::hash<int>{}(i));
		}

		auto t0 = std::chrono::high_resolution_clock::now();

		while (true) {
			render(img, pathTrace);
			img.sample_count++;

			float seconds = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t0).count() / 1000000.0);
			printf("Accumulated %d samples per pixel over %f seconds\n\n", SAMPLES_PER_PIXEL * img.sample_count, seconds);

			if (should_reset.load()) {
				should_reset = false;
				break;
			}
		}
	}
}

void process_image(sf::Image &image, std::array<std::array<Vec3, WIDTH>, HEIGHT> &img) {
	toneMap(img);
	// Load image
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			image.setPixel(px, py, img[py][px].tosRGB());
#ifdef RT_DEBUG
			Vec3 c = img[py][px];
			if (std::isnan(c.x) || std::isnan(c.y) || std::isnan(c.z)) {
				image.setPixel(px, py, sf::Color(255.0f, 0.0f, 255.0f);
		}
#endif
		}
	}
}

void gui_thread(Img &img) {
	
	sf::Image image;
	image.create(WIDTH, HEIGHT);

	//process_image(image, img);

	sf::ContextSettings contextSettings;
	contextSettings.sRgbCapable = true;

	// Init texture
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer", sf::Style::Default, contextSettings);
	sf::Texture texture;
	sf::Sprite sprite;
	texture.setSrgb(true);

	std::ofstream ofile;
	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			std::this_thread::sleep_for(100ms);
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				switch (event.key.code) {
				//case sf::Keyboard::S:
				//	printf("Saving render...\n");
				//	image.saveToFile(".\\renders\\render" + getDateTime() + ".png"); // TODO save the file as srgb
				//	ofile.open(".\\renders\\renderData" + getDateTime() + ".txt");
				//	for (int py = 0; py < HEIGHT; py++) {
				//		for (int px = 0; px < WIDTH; px++) {
				//			Vec3 v = img.pixels[py][px];
				//			ofile << (int)v.x << ',' << (int)v.y << ',' << (int)v.z << ' ';
				//		}
				//		ofile << endl;
				//	}
				//	ofile.close();

				//	printf("Render saved\n");
				//	break;
				case sf::Keyboard::Escape:
					window.close();
					break;
				default:
					break;
				}
				break;
			case sf::Event::MouseButtonPressed:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					printf("Mouse clicked at (%d, %d)\n", event.mouseButton.x, event.mouseButton.y);
				}
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
	
	auto t1 = std::chrono::high_resolution_clock::now();
	buildScene(ACTIVE_SCENE);
	float seconds = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count() / 1000000.0);
	printf("Scene built in %f seconds\n", seconds);

	Img img = Img();

	should_reset = false;

	// Render scene
	std::thread render_thread = std::thread(render_loop, std::ref(img));
	gui_thread(img);
	exit(0);
	return 0;
}