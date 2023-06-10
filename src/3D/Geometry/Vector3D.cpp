#include <math.h>
#include <ostream>

#include "Vector3D.h"

Scalar Vector3D::magSqr() const { return x * x + y * y + z * z; }

Scalar Vector3D::mag() const { return std::sqrt(x * x + y * y + z * z); }

Vector3D Vector3D::unit() const { return *this / mag(); }

Vector3D &Vector3D::operator+=(const Vector3D &rhs) {
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

Vector3D &Vector3D::operator-=(const Vector3D &rhs) {
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}

Vector3D &Vector3D::operator*=(Scalar rhs) {
  x *= rhs;
  y *= rhs;
  z *= rhs;
  return *this;
}

Vector3D &Vector3D::operator/=(Scalar rhs) {
  x /= rhs;
  y /= rhs;
  z /= rhs;
  return *this;
}

//- Functions

std::ostream &operator<<(std::ostream &os, const Vector3D &vec) {
  os << "(" << vec.x << "," << vec.y << "," << vec.z << ")";
  return os;
}

Scalar dot(const Vector3D &u, const Vector3D &v) {
  return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vector3D cross(const Vector3D &u, const Vector3D &v) {
  return Vector3D(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                  u.x * v.y - u.y * v.x);
}

Vector3D operator+(Vector3D u, const Vector3D &v) {
  u += v;
  return u;
}

Vector3D operator-(Vector3D u, const Vector3D &v) {
  u -= v;
  return u;
}

Vector3D operator*(Vector3D lhs, Scalar rhs) {
  lhs *= rhs;
  return lhs;
}

Vector3D operator*(Scalar lhs, Vector3D rhs) {
  rhs *= lhs;
  return rhs;
}

Vector3D operator/(Vector3D lhs, Scalar rhs) {
  lhs /= rhs;
  return lhs;
}
