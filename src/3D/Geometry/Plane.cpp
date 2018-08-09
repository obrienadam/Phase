#include "Plane.h"

Plane::Plane(const Point3D &pt, Vector3D n)
{
    n = n.unit();
    Scalar d = dot(n, pt);
    _coeffs = {n.x, n.y, n.z, d};
}

Plane::Plane(const Point3D &x1, const Point3D &x2, const Point3D &x3)
{
    Vector3D n = cross(x2 - x1, x3 - x1).unit();
    Scalar d = dot(n, x1);
    _coeffs = {n.x, n.y, n.z, d};
}

Point3D Plane::nearestPoint(const Point3D &pt) const
{
    Scalar t = (_coeffs[3] - dot(n(), pt)) / n().magSqr();
    return pt + t * n();
}
