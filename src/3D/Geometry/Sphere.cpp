#include "Sphere.h"

Sphere::Sphere(const Point3D &centroid, Scalar radius)
    :
      Shape3D(centroid),
      radius_(radius)
{
    volume_ = 4. * M_PI * radius_ * radius_ * radius_ / 3.;
}

bool Sphere::isInside(const Point3D &pt) const
{
    return (pt - centroid_).magSqr() < radius_ * radius_;
}
