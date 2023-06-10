#include "Ray2D.h"

Ray2D::Ray2D(const Point2D &x0, const Vector2D &r) : x0_(x0), r_(r.unitVec()) {}

Point2D Ray2D::operator()(Scalar t) const { return x0_ + r_ * t; }

bool Ray2D::isBounded(const Point2D &pt) const {
  return dot(pt - x0_, r_) >= 0.;
}

Ray2D Ray2D::rotate(Scalar theta) const { return Ray2D(x0_, r_.rotate(theta)); }