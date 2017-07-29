#ifndef TENSOR_2D_H
#define TENSOR_2D_H

#include "Types.h"
#include "Vector2D.h"

class Tensor2D
{
public:

    Tensor2D(Scalar xx = 0., Scalar xy = 0., Scalar yx = 0., Scalar yy = 0.) : xx(xx), xy(xy), yx(yx), yy(yy)
    {}

    Scalar xx, xy, yx, yy;

    Tensor2D &operator+=(const Tensor2D &rhs);

    Tensor2D &operator-=(const Tensor2D &rhs);

    Tensor2D &operator*=(Scalar a);

    Tensor2D &operator/=(Scalar a);

};

std::ostream &operator<<(std::ostream &os, const Tensor2D &tau);

Vector2D dot(const Tensor2D &tau, const Vector2D &u);

Tensor2D dot(const Tensor2D &tau, const Tensor2D &sigma);

Tensor2D outer(const Vector2D &u, const Vector2D &v);

Tensor2D operator+(Tensor2D lhs, const Tensor2D &rhs);

Tensor2D operator-(Tensor2D lhs, const Tensor2D &rhs);

Tensor2D operator*(Tensor2D lhs, Scalar a);

Tensor2D operator*(Scalar a, Tensor2D rhs);

Tensor2D operator/(Tensor2D lhs, Scalar a);

#endif
