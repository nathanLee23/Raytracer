#include <math.h>
#include "Matrix3.h"

Vec3 operator *(const Matrix3& m, const Vec3 v) {
	Vec3 mv = Vec3();
	mv.x += m.matrix[0][0] * v.x + m.matrix[0][1] * v.y + m.matrix[0][2] * v.z;
	mv.y += m.matrix[1][0] * v.x + m.matrix[1][1] * v.y + m.matrix[1][2] * v.z;
	mv.z += m.matrix[2][0] * v.x + m.matrix[2][1] * v.y + m.matrix[2][2] * v.z;

	return mv;
}

//Creates a matrix that rotates normalized vectors b to a
Matrix3 rotMatrixVectors(Vec3 a, Vec3 b) {
	Vec3 axis = b.cross(a);

	float c = a.dot(b);
	float s = axis.norm();

	Matrix3 rotMatrix;
	rotMatrix.matrix[0][0] = c + axis.x*axis.x*(1.0f - c);
	rotMatrix.matrix[0][1] = axis.x*axis.y*(1.0f - c) - axis.z*s;
	rotMatrix.matrix[0][2] = axis.x*axis.z*(1.0f - c) + axis.y*s;
	rotMatrix.matrix[1][0] = rotMatrix.matrix[0][1] + 2 * axis.z*s;
	rotMatrix.matrix[1][1] = c + axis.y*axis.y*(1.0f - c);
	rotMatrix.matrix[1][2] = axis.y*axis.z*(1.0f - c) - axis.x*s;
	rotMatrix.matrix[2][0] = rotMatrix.matrix[0][2] - 2 * axis.y*s;
	rotMatrix.matrix[2][1] = rotMatrix.matrix[1][2] + 2 * axis.x*s;
	rotMatrix.matrix[2][2] = c + axis.z*axis.z*(1.0f - c);
	return rotMatrix;
}

Matrix3 rotMatrixOnAxis(const Vec3 axis, float angle) {
	float c = std::cos(angle);
	float s = std::sin(angle);

	Matrix3 rotMatrix;
	rotMatrix.matrix[0][0] = c + axis.x*axis.x*(1.0f - c);
	rotMatrix.matrix[0][1] = axis.x*axis.y*(1.0f - c) - axis.z*s;
	rotMatrix.matrix[0][2] = axis.x*axis.z*(1.0f - c) + axis.y*s;
	rotMatrix.matrix[1][0] = rotMatrix.matrix[0][1] + 2 * axis.z*s;
	rotMatrix.matrix[1][1] = c + axis.y*axis.y*(1.0f - c);
	rotMatrix.matrix[1][2] = axis.y*axis.z*(1.0f - c) - axis.x*s;
	rotMatrix.matrix[2][0] = rotMatrix.matrix[0][2] - 2 * axis.y*s;
	rotMatrix.matrix[2][1] = rotMatrix.matrix[1][2] + 2 * axis.x*s;
	rotMatrix.matrix[2][2] = c + axis.z*axis.z*(1.0f - c);
	return rotMatrix;
}