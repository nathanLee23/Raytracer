#include <iostream>
#include "Vec3.h"
Vec3::Vec3() : x(0.0), y(0.0), z(0.0) {
}

Vec3::Vec3(double x) : x(x), y(x), z(x) {
}

Vec3::Vec3(double x, double y, double z) : x(x), y(y), z(z) {
}

double Vec3::dot(const Vec3 o) {
	return x * o.x + y * o.y + z * o.z;
}
Vec3 Vec3::cross(const Vec3 o) {
	return Vec3(y*o.z - z * o.y, z*o.x - x * o.z, x*o.y - y * o.x);
}
double Vec3::sqrNorm() {
	return x * x + y * y + z * z;
}
double Vec3::norm() {
	return sqrt(x * x + y * y + z * z);
}
Vec3 Vec3::normalized() {
	double n = norm();
	return Vec3(x / n, y / n, z / n);
}


Vec3 operator+(const Vec3 a, const Vec3 b) {
	return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vec3 operator*(const Vec3 a, double b) {
	return Vec3(a.x *b, a.y *b, a.z *b);
}
Vec3 operator*(double b, const Vec3 a) {
	return a * b;
}
Vec3 operator/(const Vec3 a, double b) {
	return a * (1 / b);
}
Vec3 operator-(const Vec3 a) {
	return Vec3(-a.x, -a.y, -a.z);
}
Vec3 operator-(const Vec3 a, const Vec3 b) {
	return a + (-b);
}
std::ostream& operator<<(std::ostream& os, const Vec3 a) {
	os << "(" << a.x << ", " << a.y << ", " << a.z << ")";
	return os;
}