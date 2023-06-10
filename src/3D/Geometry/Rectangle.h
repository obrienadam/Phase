#ifndef PHASE_RECTANGLE_H
#define PHASE_RECTANGLE_H

#include "Shape2D.h"

class Rectangle : public Shape2D {
public:
  Rectangle(const Point3D &centroid = Point3D(0., 0., 0.),
            const Vector3D &norm = Vector3D(0., 0., 0.));

  Point3D project(const Point3D &pt) const override;

  const Vector3D &norm() const { return norm_; }

private:
  Vector3D norm_;
};

#endif
