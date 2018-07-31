#ifndef PHASE_VECTOR_3D_H
#define PHASE_VECTOR_3D_H

#include "Types/Types.h"

class Vector3D
{
public:

    Vector3D(Scalar x = 0., Scalar y = 0., Scalar z = 0.) : x(x), y(y), z(z) {}

    Scalar magSqr() const;

    Scalar mag() const;

    Vector3D unit() const;

    Vector3D &operator+=(const Vector3D &rhs);

    Vector3D &operator-=(const Vector3D &rhs);

    Vector3D &operator*=(Scalar rhs);

    Vector3D &operator/=(Scalar rhs);

    Scalar x, y, z;
};

std::ostream &operator<<(std::ostream &os, const Vector3D &vec);

Scalar dot(const Vector3D &u, const Vector3D &v);

Vector3D cross(const Vector3D &u, const Vector3D &v);

Vector3D operator+(Vector3D u, const Vector3D &v);

Vector3D operator-(Vector3D u, const Vector3D &v);

Vector3D operator*(Vector3D lhs, Scalar rhs);

Vector3D operator*(Scalar lhs, Vector3D rhs);

Vector3D operator/(Vector3D lhs, Scalar rhs);

#endif
