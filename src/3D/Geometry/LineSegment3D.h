#ifndef PHASE_LINE_SEGMENT_3D_H
#define PHASE_LINE_SEGMENT_3D_H

#include "Point3D.h"

class LineSegment3D {
public:
  LineSegment3D(const Point3D &ptA, const Point3D &ptB)
      : _ptA(ptA), _ptB(ptB) {}

  const Point3D &ptA() const { return _ptA; }

  const Point3D &ptB() const { return _ptB; }

  Vector3D r() const { return _ptB - _ptA; }

  Scalar length() const { return (_ptB - _ptA).mag(); }

  Scalar lengthSqr() const { return (_ptB - _ptA).magSqr(); }

private:
  Point3D _ptA, _ptB;
};

#endif
