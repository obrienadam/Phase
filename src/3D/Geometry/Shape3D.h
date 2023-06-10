#ifndef PHASE_SHAPE_3D_H
#define PHASE_SHAPE_3D_H

#include "Point3D.h"

class Shape3D {
public:
  Shape3D(const Point3D &centroid = Point3D(0., 0., 0.))
      : centroid_(centroid) {}

  Scalar volume() const { return volume_; }

  const Point3D &centroid() const { return centroid_; }

  //- Tests
  virtual bool isInside(const Point3D &pt) const = 0;

protected:
  Scalar volume_;

  Point3D centroid_;
};

#endif
