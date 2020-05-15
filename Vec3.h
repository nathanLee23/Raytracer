#pragma once
#include "SFML/Graphics/Color.hpp"

#include <iostream>
class Vec3 {
public:
	double x, y, z;
	Vec3();
	Vec3(double x);
	Vec3(double x, double y, double z);

	double dot(const Vec3 o);
	Vec3 cross(const Vec3 o);
	double sqrNorm();
	double norm();
	Vec3 normalized();
	Vec3 proj(const Vec3 a); // Projects vector a onto this
	sf::Color toColor();
};


Vec3 operator+(const Vec3 a, const Vec3 b);
Vec3 operator*(const Vec3 a, double b);
Vec3 operator*(double b, const Vec3 a);
Vec3 operator/(const Vec3 a, double b);
Vec3 operator-(const Vec3 a);
Vec3 operator-(const Vec3 a, const Vec3 b);
std::ostream& operator<<(std::ostream& os, const Vec3 a);