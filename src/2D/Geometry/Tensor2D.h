#ifndef TENSOR_2D_H
#define TENSOR_2D_H

#include "Vector2D.h"

class Tensor2D
{
public:

    static Tensor2D eye()
    { return Tensor2D(1., 0., 0., 1.); }

    Tensor2D(Scalar xx = 0., Scalar xy = 0., Scalar yx = 0., Scalar yy = 0.) : xx(xx), xy(xy), yx(yx), yy(yy)
    {}

    Scalar xx, xy, yx, yy;

    Tensor2D transpose() const
    { return Tensor2D(xx, yx, xy, yy); }

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

Tensor2D operator-(const Tensor2D &tau);

Tensor2D operator*(Tensor2D lhs, Scalar a);

Tensor2D operator*(Scalar a, Tensor2D rhs);

Tensor2D operator/(Tensor2D lhs, Scalar a);

#endif
