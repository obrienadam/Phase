#include "Tensor2D.h"

//- Public methods

Tensor2D &Tensor2D::operator+=(const Tensor2D &rhs) {
  xx += rhs.xx;
  xy += rhs.xy;
  yx += rhs.yx;
  yy += rhs.yy;

  return *this;
}

Tensor2D &Tensor2D::operator-=(const Tensor2D &rhs) {
  xx -= rhs.xx;
  xy -= rhs.xy;
  yx -= rhs.yx;
  yy -= rhs.yy;

  return *this;
}

Tensor2D &Tensor2D::operator*=(Scalar a) {
  xx *= a;
  xy *= a;
  yx *= a;
  yy *= a;

  return *this;
}

Tensor2D &Tensor2D::operator/=(Scalar a) {
  xx /= a;
  xy /= a;
  yx /= a;
  yy /= a;

  return *this;
}

std::ostream &operator<<(std::ostream &os, const Tensor2D &tau) {
  os << "[xx=" << tau.xx << " xy=" << tau.xy << "; yx=" << tau.yx
     << " yy=" << tau.yy << "]";
  return os;
}

//- External functions
Vector2D dot(const Tensor2D &tau, const Vector2D &u) {
  return Vector2D(tau.xx * u.x + tau.xy * u.y, tau.yx * u.x + tau.yy * u.y);
}

Tensor2D dot(const Tensor2D &tau, const Tensor2D &sigma) {
  return Tensor2D(tau.xx * sigma.xx + tau.xy * sigma.yx,
                  tau.xx * sigma.xy + tau.xy * sigma.yy,
                  tau.yx * sigma.xx + tau.yy * sigma.yx,
                  tau.yx * sigma.xy + tau.yy * sigma.yy);
}

Tensor2D outer(const Vector2D &u, const Vector2D &v) {
  return Tensor2D(u.x * v.x, u.x * v.y, u.y * v.x, u.y * v.y);
}

Tensor2D operator+(Tensor2D lhs, const Tensor2D &rhs) { return lhs += rhs; }

Tensor2D operator-(Tensor2D lhs, const Tensor2D &rhs) { return lhs -= rhs; }

Tensor2D operator-(const Tensor2D &tau) {
  return Tensor2D(-tau.xx, -tau.xy, -tau.yx, -tau.yy);
}

Tensor2D operator*(Tensor2D lhs, Scalar a) { return lhs *= a; }

Tensor2D operator*(Scalar a, Tensor2D rhs) { return rhs *= a; }

Tensor2D operator/(Tensor2D lhs, Scalar a) { return lhs /= a; }
