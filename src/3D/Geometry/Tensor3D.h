#ifndef PHASE_TENSOR_3D_H
#define PHASE_TENSOR_3D_H

#include "Vector3D.h"

class Tensor3D {
public:
  Tensor3D(Scalar xx = 0., Scalar xy = 0., Scalar xz = 0., Scalar yx = 0.,
           Scalar yy = 0., Scalar yz = 0., Scalar zx = 0., Scalar zy = 0.,
           Scalar zz = 0.)
      : xx(xx), xy(xy), xz(xz), yx(yx), yy(yy), yz(yz), zx(zx), zy(zy), zz(zz) {
  }

  Scalar xx, xy, xz, yx, yy, yz, zx, zy, zz;

  //- Operators
  Tensor3D &operator+=(const Tensor3D &tau);

  Tensor3D &operator-=(const Tensor3D &tau);

  Tensor3D &operator*=(Scalar a);

  Tensor3D &operator/=(Scalar a);

  Vector3D trace() const { return Vector3D(xx, yy, zz); }

  Tensor3D transpose() const;
};

Vector3D dot(const Tensor3D &tau, const Vector3D &u);

Tensor3D outer(const Vector3D &u, const Vector3D &v);

#endif
