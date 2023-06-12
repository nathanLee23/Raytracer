// Raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <thread>
#include <functional>
#include <string.h>
#include <atomic>

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <random> // TODO Replace with QMC generators

#include <chrono>
#include <ctime>  

#include <omp.h>
#include <SFML/Graphics.hpp>
#include <embree3/rtcore.h>

#define TINYEXR_IMPLEMENTATION 
#include "tinyexr.h"

#include "Vec3.h"
#include "Matrix3.h"
#include "Material.h"

#include "Obj.h"
#include "globals.h"
#include "Camera.h"
#include "Scene.h"
#include "Integrator.h"

#define FOVX 60.0f
#define ACTIVE_SCENE 10

#define AMBIENT 0.3f
#define SAMPLES_PER_PIXEL 1
#define THREADS 8

using namespace std::chrono_literals;

std::atomic<bool> should_reset;

class Img {
public:
	std::vector<Vec3> pixels;

	unsigned int sample_count;
	Img() {
		pixels.resize(WIDTH*HEIGHT);
		clear();
	}
	
	static uint32_t idx(int px, int py) {
		return px + py * WIDTH;
	}

	void clear() {
		for (int py = 0; py < HEIGHT; py++) {
			for (int px = 0; px < WIDTH; px++) {
				pixels[idx(px, py)] = Vec3();
			}
		}
		sample_count = 0u;
	}	
	
	void set(int px, int py, Vec3 v) {
		pixels[idx(px, py)] = v;
	}

	Vec3 get(int px, int py) const {
		return pixels[idx(px, py)];
	}

	void update(int px, int py, Vec3 color) {
		// TODO use Kahan summation?
		set(px, py, (get(px, py) * static_cast<float>(sample_count) + color) / static_cast<float>(sample_count + 1));
	}

	bool saveEXR(const char* outfilename) {

		EXRHeader header;
		InitEXRHeader(&header);

		EXRImage image;
		InitEXRImage(&image);

		image.num_channels = 3;

		std::vector<float> images[3];
		images[0].resize(WIDTH * HEIGHT);
		images[1].resize(WIDTH * HEIGHT);
		images[2].resize(WIDTH * HEIGHT);

		// Split RGBRGBRGB... into R, G and B layer
		for (int i = 0; i < HEIGHT; i++) {
			for (int j = 0; j < WIDTH; j++) {
				images[0][j * HEIGHT + i] = get(i, j).x;
				images[1][j * HEIGHT + i] = get(i, j).y;
				images[2][j * HEIGHT + i] = get(i, j).z;
			}
		}

		float* image_ptr[3];
		image_ptr[0] = &(images[2].at(0)); // B
		image_ptr[1] = &(images[1].at(0)); // G
		image_ptr[2] = &(images[0].at(0)); // R

		image.images = (unsigned char**)image_ptr;
		image.width = WIDTH;
		image.height = HEIGHT;

		header.num_channels = 3;
		header.channels = (EXRChannelInfo*)malloc(sizeof(EXRChannelInfo) * header.num_channels);
		// Must be (A)BGR order, since most of EXR viewers expect this channel order.
		strncpy(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
		strncpy(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
		strncpy(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';

		header.pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
		header.requested_pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
		for (int i = 0; i < header.num_channels; i++) {
			header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
			header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
		}

		const char* err = NULL; // or nullptr in C++11 or later.
		int ret = SaveEXRImageToFile(&image, &header, outfilename, &err);
		if (ret != TINYEXR_SUCCESS) {
			fprintf(stderr, "Save EXR err: %s\n", err);
			FreeEXRErrorMessage(err); // free's buffer for an error message
			return ret;
		}
		printf("Saved exr file. [ %s ] \n", outfilename);

		free(header.channels);
		free(header.pixel_types);
		free(header.requested_pixel_types);

	}
};
auto camera = PerspectiveCamera(50.0f);
//auto camera = OrthogonalCamera(1.2f);
//auto camera = ThinLensCamera();
Scene scene = Scene();

void buildScene(int i) {
	Material mirrorMat = { Vec3(1.0f), 0.0f, Surface(reflective) };
	Material diffuseMat = { Vec3(0.93f, 0.93f, 0.93f), 0.0f, Surface(diffuse) };
	Material varnishMat = { Vec3(0.73f, 0.73f, 0.73f), 0.0f, Surface(varnish) };
	Material redMat = { Vec3(0.55f, 0.09f, 0.09f), 0.0f, Surface(diffuse) };
	Material greenMat = { Vec3(0.16f, 0.55f, 0.15f), 0.0f, Surface(diffuse) };
	Material specularMat = { Vec3(1.0f), 0.0f, Surface(specular) };

	// Using specular light sauces creates a lot of noise
	Material lightMat = { Vec3(1.0f), 16.0f, Surface(specular) };
	Material ovenMat = { Vec3(0.5f), 0.5f, Surface(diffuse) };

	scene.addMesh(scene.load_mesh("geometry/CornellBox-Original.obj"));
	scene.meshes[0]->materials[0] = redMat; // left wall
	scene.meshes[0]->materials[1] = greenMat; // right wall
	scene.meshes[0]->materials[2] = diffuseMat; // floor
	scene.meshes[0]->materials[3] = diffuseMat;
	scene.meshes[0]->materials[4] = diffuseMat;
	scene.meshes[0]->materials[5] = diffuseMat; // right box
	scene.meshes[0]->materials[6] = diffuseMat; //left box
	//scene.meshes[0]->materials[6] = mirrorMat; //left box
	scene.meshes[0]->materials[7].emission = 10.0f; // light

#ifdef SPHERES
	scene.spheres.push_back(Sphere(Vec3(-0.5f,0.302f,0.55f), 0.3f));
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

Vec3 pos;
template<typename I>
void render(Img &img, I&& integrator) {
	// TODO Use tiles, because we assign threads adjacent pixels we are suffering from false sharing of the image buffer

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
					ray.o += pos;
					colorVec += integrator(scene, ray, gen);
				}
				img.set(px, py, colorVec/SAMPLES_PER_PIXEL);
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
					ray.o += pos;
					colorVec += integrator(scene, ray, gen);
				}
				img.update(px, py, colorVec/SAMPLES_PER_PIXEL);
			}
		}

	}
		
	// Performance metrics
	float seconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - t1).count();
	printf("Took %f seconds\n", seconds);
	printf("Cast %d rays\n", SAMPLES_PER_PIXEL * WIDTH*HEIGHT);
	printf("%f samples per second \n", SAMPLES_PER_PIXEL * WIDTH * HEIGHT / seconds);
	printf("%f ms per sample \n", seconds * 1'000'000.0f / (SAMPLES_PER_PIXEL * WIDTH * HEIGHT));
}

void toneMap(Img &img) {
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
	// TODO Schedule work adaptively like how blender/cycles does it: https://developer.blender.org/diffusion/C/browse/master/src/integrator/render_scheduler.cpp$809

	while (true) {
		for (int i = 0; i < THREADS; i++) {
			gens[i] = Pcg(std::hash<int>{}(i));
		}

		img.clear();
		auto t0 = std::chrono::high_resolution_clock::now();

		while (true) {
			render(img, pathTrace);
			img.sample_count++;

			float seconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - t0).count();
			printf("Accumulated %d samples per pixel over %f seconds\n\n", SAMPLES_PER_PIXEL * img.sample_count, seconds);

			if (should_reset.load()) {
				should_reset = false;
				break;
			}
		}
	}
}

void process_image(sf::Image &image, Img &img) {
	toneMap(img);
	// Load image
	for (int py = 0; py < HEIGHT; py++) {
		for (int px = 0; px < WIDTH; px++) {
			image.setPixel(px, py, img.get(px, py).tosRGB());
#ifdef RT_DEBUG
			Vec3 c = img.get(px, py);
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
	sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Ray tracer", sf::Style::Titlebar | sf::Style::Close, contextSettings);
	sf::Texture texture;
	sf::Sprite sprite;
	texture.setSrgb(true);

	std::ofstream ofile;
	std::string save_loc;
	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				switch (event.key.code) {
				case sf::Keyboard::Space:
					save_loc = ".\\renders\\render" + getDateTime() + ".exr";
					img.saveEXR(save_loc.data());
					break;
				case sf::Keyboard::Escape:
					window.close();
					break;
				case sf::Keyboard::R:
					should_reset = true;
					break;
				default:
					break;
				}
				break;
			case sf::Event::MouseButtonPressed:
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					Vec3 col = img.get(event.mouseButton.x, event.mouseButton.y);
					printf("(%f, %f, %f) at (%d, %d)\n", col.x, col.y, col.z, event.mouseButton.x, event.mouseButton.y);
				}
			default:
				break;
			}
		}

		// TODO Camera rotation by precomputing the projection plane!!
		if (window.hasFocus())
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
			{
				pos += Vec3(1.0, 0.0, 0.0) * 0.01;
				should_reset = true;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
			{
				pos += Vec3(-1.0, 0.0, 0.0) * 0.01;
				should_reset = true;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
			{
				pos += Vec3(0.0, 0.0, -1.0) * 0.01;
				should_reset = true;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
			{
				pos += Vec3(0.0, 0.0, 1.0) * 0.01;
				should_reset = true;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift))
			{
				pos += Vec3(0.0, 1.0, 0.0) * 0.01;
				should_reset = true;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::LControl))
			{
				pos += Vec3(0.0, -1.0, 0.0) * 0.01;
				should_reset = true;
			}
		}

		window.clear();
		window.draw(sprite);
		window.display();

		process_image(image, img);
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
	float seconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - t1).count();
	printf("Scene built in %f seconds\n", seconds);

	Img img = Img();

	should_reset = false;

	// Render scene
	std::thread render_thread = std::thread(render_loop, std::ref(img));
	gui_thread(img);
	exit(0);
	return 0;
}