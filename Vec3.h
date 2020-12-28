#pragma once
#include <iostream>

#include "SFML/Graphics/Color.hpp"
class Vec3 {
public:
	float x, y, z;
	Vec3();
	Vec3(float x);
	Vec3(float x, float y, float z);

	float dot(const Vec3 o);
	Vec3 cross(const Vec3 o);
	float sqrNorm();
	float norm();
	float max();
	Vec3 normalized();
	Vec3 proj(const Vec3 a); // Projects vector a onto this
	sf::Color tosRGB();

	float operator[] (int i);
};


Vec3 operator+(const Vec3 a, const Vec3 b);
Vec3 operator*(const Vec3 a, float b);
Vec3 operator*(float b, const Vec3 a);
Vec3 operator*(const Vec3 a, const Vec3 b);
Vec3 operator/(const Vec3 a, float b);
Vec3 operator-(const Vec3 a);
Vec3 operator-(const Vec3 a, const Vec3 b);
std::ostream& operator<<(std::ostream& os, const Vec3 a);