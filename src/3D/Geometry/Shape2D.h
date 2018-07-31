#ifndef PHASE_SHAPE_2D_H
#define PHASE_SHAPE_2D_H

#include "Point3D.h"

class Shape2D
{
public:

    Shape2D(const Point3D &centroid = Point3D(0., 0., 0.)) : centroid_(centroid) {}

    Scalar area() const
    { return area_; }

    const Point3D &centroid() const
    { return centroid_; }

    virtual Point3D project(const Point3D &pt) const = 0;

protected:

    Scalar area_;

    Point3D centroid_;

};

#endif
