#include <iostream>
#include <algorithm>

#include "Vec3.h"

Vec3::Vec3() : x(0.0f), y(0.0f), z(0.0) {
}

Vec3::Vec3(float x) : x(x), y(x), z(x) {
}

Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {
}

float Vec3::dot(const Vec3 o) const {
	return x * o.x + y * o.y + z * o.z;
}
Vec3 Vec3::cross(const Vec3 o) const {
	return Vec3(y*o.z - z * o.y, z*o.x - x * o.z, x*o.y - y * o.x);
}
float Vec3::sqrNorm() const {
	return x * x + y * y + z * z;
}
float Vec3::norm() const {
	return std::sqrt(x * x + y * y + z * z);
}
Vec3 Vec3::normalized() const {
	float n = norm();
	return Vec3(x / n, y / n, z / n);
}
Vec3 Vec3::proj(const Vec3 a) const {
	return dot(a) / sqrNorm() * (*this);
}

float Vec3::max() const {
	return std::max(x, std::max(y, z));
}

float correctGamma(float x) {
	if (x <= 0.0031308f) {
		return 12.92f * x;
	} else {
		return 1.055f * std::powf(x, 1 / 2.4) - 0.055;
	}
}

sf::Color Vec3::tosRGB() const {
	return sf::Color(
		255.0f*correctGamma(std::min(1.0f, x)),
		255.0f*correctGamma(std::min(1.0f, y)),
		255.0f*correctGamma(std::min(1.0f, z))
	);
}

sf::Color Vec3::toColor() const {
	return sf::Color(
		255.0f * std::min(1.0f, x),
		255.0f * std::min(1.0f, y),
		255.0f * std::min(1.0f, z)
	);
}

Vec3 operator+(const Vec3& a, const Vec3& b) {
	return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vec3 operator*(const Vec3& a, const float b) {
	return Vec3(a.x *b, a.y *b, a.z *b);
}
Vec3 operator*(float b, const Vec3& a) {
	return a * b;
}
Vec3 operator*(const Vec3& a, const Vec3& b) {
	return Vec3(a.x*b.x, a.y*b.y, a.z*b.z);
}
Vec3 operator/(const Vec3& a, const Vec3& b) {
	return Vec3(a.x / b.x, a.y / b.y, a.z / b.z);
}
Vec3 operator/(const Vec3& a, const float b) {
	return a * (1.0f / b);
}
Vec3 operator/(const float a, const Vec3& b) {
	return Vec3(a / b.x, a / b.y, a / b.z);
}
Vec3 operator-(const Vec3& a) {
	return Vec3(-a.x, -a.y, -a.z);
}
Vec3 operator-(const Vec3& a, const Vec3& b) {
	return a + (-b);
}

Vec3& operator+=(Vec3& a, const Vec3& b) {
	a = a + b;
	return a;
}

Vec3& operator-=(Vec3& a, const Vec3& b) {
	a = a - b;
	return a;
}

Vec3& operator*=(Vec3& a, const Vec3& b) {
	a = a * b;
	return a;
}

Vec3& operator*=(Vec3& a, const float b) {
	a = a * b;
	return a;
}

Vec3& operator/=(Vec3& a, const float b) {
	a = a / b;
	return a;
}


std::ostream& operator<<(std::ostream& os, const Vec3 a) {
	os << "(" << a.x << ", " << a.y << ", " << a.z << ")";
	return os;
}