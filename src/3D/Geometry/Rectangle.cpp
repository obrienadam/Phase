#include "Rectangle.h"

Rectangle::Rectangle(const Point3D &centroid, const Vector3D &norm)
    : norm_(norm) {
  centroid_ = centroid;
  area_ = norm.mag();
}

Point3D Rectangle::project(const Point3D &pt) const {
  Vector3D r = pt - centroid_;
  return r - dot(r, norm_.unit()) + centroid_;
}
