#pragma once
#include <iostream>

#include "SFML/Graphics/Color.hpp"
class Vec3 {
public:
	float x, y, z;
	Vec3();
	Vec3(float x);
	Vec3(float x, float y, float z);

	float dot(const Vec3 o) const;
	Vec3 cross(const Vec3 o) const;
	float sqrNorm() const;
	float norm() const;
	float max() const;
	Vec3 normalized() const;
	Vec3 proj(const Vec3 a) const; // Projects vector a onto this
	sf::Color tosRGB() const;
	sf::Color toColor() const;
};


Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator*(const Vec3& a, const float b);
Vec3 operator*(float b, const Vec3& a);
Vec3 operator*(const Vec3& a, const Vec3& b);
Vec3 operator/(const Vec3& a, const Vec3& b);
Vec3 operator/(const Vec3& a, const float b);
Vec3 operator-(const Vec3& a);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3& operator+=(Vec3& a, const Vec3& b);
Vec3& operator-=(Vec3& a, const Vec3& b);
Vec3& operator*=(Vec3& a, const Vec3& b);
Vec3& operator*=(Vec3& a, const float b);
Vec3& operator/=(Vec3& a, const float b);
std::ostream& operator<<(std::ostream& os, const Vec3 a);