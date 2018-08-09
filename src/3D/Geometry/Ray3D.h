#ifndef PHASE_RAY_3D_H
#define PHASE_RAY_3D_H

#include "Point3D.h"

class Ray3D
{
public:

    Ray3D(const Point3D &x0 = Point3D(0., 0., 0.), const Vector3D &r = Vector3D(0., 0., 0.)) : _x0(x0), _r(r) {}

    const Point3D &x0() const
    { return _x0; }

    const Vector3D &r() const
    { return _r; }

    Point3D operator()(Scalar t) const
    { return _x0 + t * _r; }

private:

    Point3D _x0;

    Vector3D _r;
};

#endif
