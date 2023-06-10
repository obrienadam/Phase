#ifndef PHASE_RECTANGULAR_PRISM_H
#define PHASE_RECTANGULAR_PRISM_H

#include "Shape3D.h"

class RectangularPrism : public Shape3D {
public:
  RectangularPrism(const Point3D &minVertex = Point3D(0., 0., 0.),
                   const Point3D &maxVertex = Point3D(0., 0., 0.));

  bool isInside(const Point3D &pt) const override;

private:
  Point3D minVertex_, maxVertex_;
};

#endif
