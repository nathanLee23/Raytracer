#pragma once

#include "Vec3.h"

class Matrix3 {
public:
	double matrix[3][3];
};

Vec3 operator *(const Matrix3& m, const Vec3 v);
Matrix3 rotMatrixVectors(Vec3 a, Vec3 b);
Matrix3 rotMatrixOnAxis(const Vec3 axis, double angle);