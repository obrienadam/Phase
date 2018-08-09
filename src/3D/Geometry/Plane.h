#ifndef PHASE_PLANE_H
#define PHASE_PLANE_H

#include "Point3D.h"

class Plane
{
public:

    Plane() {}

    Plane(const Point3D &pt, Vector3D n);

    Plane(const Point3D &x1, const Point3D &x2, const Point3D &x3);

    const std::array<Scalar, 4> &coeffs() const
    { return _coeffs; }

    Vector3D n() const
    { return Vector3D(_coeffs[0], _coeffs[1], _coeffs[2]); }

    Point3D nearestPoint(const Point3D &pt) const;

private:

    std::array<Scalar, 4> _coeffs;

};

#endif
