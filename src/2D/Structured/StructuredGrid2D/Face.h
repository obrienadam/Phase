#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "Cell.h"

class Face {
public:
  Face(const StructuredGrid2D &grid, Coordinates::Index coord, Label i,
       Label j);

  const LineSegment2D &shape() const { return _shape; }

  Point2D centroid() const { return _shape.center(); }

  Label id() const { return _lid; }

  Label lid() const { return _lid; }

  const Cell &lCell() const { return *_lCell; }

  const Cell &rCell() const { return *_rCell; }

  bool isBoundaryFace() const { return (bool)_rCell; }

  const Vector2D &norm() const { return _norm; }

  Vector2D sf(const Point2D &pt) const {
    return dot(_norm, _shape.center() - pt) > 0. ? _norm : -_norm;
  }

protected:
  LineSegment2D _shape;

  Vector2D _norm;

  Label _i, _j, _lid;

  const StructuredGrid2D &_grid;

  const Cell *_lCell, *_rCell;
};

#endif
