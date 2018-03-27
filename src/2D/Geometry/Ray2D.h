#ifndef RAY_2D_H
#define RAY_2D_H

#include "Point2D.h"

class Ray2D
{
public:

    Ray2D() {}

    Ray2D(const Point2D& x0, const Vector2D& r);

    const Point2D& x0() const
    { return x0_; }

    const Vector2D& r() const
    { return r_; }

    Point2D operator()(Scalar t) const;

    bool isBounded(const Point2D& pt) const;

    Ray2D rotate(Scalar theta) const;

private:

    Point2D x0_;
    Vector2D r_;
};

#endif