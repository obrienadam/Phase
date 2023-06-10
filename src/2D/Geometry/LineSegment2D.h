#ifndef PHASE_LINE_SEGMENT_2D_H
#define PHASE_LINE_SEGMENT_2D_H

#include "Point2D.h"

class Shape2D;

class LineSegment2D {
public:
  LineSegment2D(const Point2D &ptA = Point2D(), const Point2D &ptB = Point2D())
      : ptA_(ptA), ptB_(ptB) {}

  const Point2D &ptA() const { return ptA_; }

  const Point2D &ptB() const { return ptB_; }

  Scalar lengthSqr() const { return (ptB_ - ptA_).magSqr(); }

  Scalar length() const { return (ptB_ - ptA_).mag(); }

  Scalar dx() const { return ptB_.x - ptA_.x; }

  Scalar dy() const { return ptB_.y - ptA_.y; }

  Vector2D rVec() const { return ptB_ - ptA_; }

  Point2D center() const { return (ptA_ + ptB_) / 2.; }

  Vector2D norm() const { return (ptB_ - ptA_).normalVec(); }

  bool isBounded(const Point2D &pt) const;

private:
  Point2D ptA_, ptB_;
};

#endif
