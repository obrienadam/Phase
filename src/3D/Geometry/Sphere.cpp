#include "Sphere.h"

Sphere::Sphere(const Point3D &centroid, Scalar radius)
    :
      Shape3D(centroid),
      _radius(radius)
{
    volume_ = 4. * M_PI * std::pow(_radius, 3) / 3.;
}

bool Sphere::isInside(const Point3D &pt) const
{
    return (pt - centroid_).magSqr() < std::pow(_radius, 2);
}
