#ifndef PHASE_SPHERE_H
#define PHASE_SPHERE_H

#include "Shape3D.h"

class Sphere: public Shape3D
{
public:

    Sphere(const Point3D &centroid = Point3D(0., 0., 0.), Scalar radius = 0.);

    //- Tests
    bool isInside(const Point3D &pt) const override;

private:

    Scalar radius_;
};

#endif
