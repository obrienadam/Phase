#include "LineSegment2D.h"

bool LineSegment2D::isBounded(const Point2D &pt) const {
  Scalar t = dot(pt - ptA_, ptB_ - ptA_);
  return t >= 0. && t <= lengthSqr();
}
