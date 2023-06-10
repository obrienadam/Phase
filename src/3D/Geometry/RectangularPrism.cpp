#include "RectangularPrism.h"

RectangularPrism::RectangularPrism(const Point3D &minVertex,
                                   const Point3D &maxVertex)
    : minVertex_(minVertex), maxVertex_(maxVertex) {
  volume_ = (maxVertex_.x - minVertex_.x) * (maxVertex_.y - minVertex_.y) *
            (maxVertex_.z - minVertex_.z);
  centroid_ = (minVertex_ + maxVertex_) / 2.;
}

bool RectangularPrism::isInside(const Point3D &pt) const {
  return (pt.x > minVertex_.x && pt.x < maxVertex_.x) &&
         (pt.y > minVertex_.y && pt.y < maxVertex_.y) &&
         (pt.z > minVertex_.z && pt.z < maxVertex_.z);
}
