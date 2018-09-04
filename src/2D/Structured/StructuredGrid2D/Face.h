#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "Cell.h"

class Face
{
public:

    Face(const StructuredGrid2D &grid, Coordinate coord, Label i, Label j);

    const LineSegment2D &shape() const
    { return _shape; }

    Point2D centroid() const
    { return _shape.center(); }

    Label lid() const
    { return _lid; }

    const Cell &lCell() const
    { return *_lCell; }

    const Cell &rCell() const
    { return *_rCell; }

    bool isBoundaryFace() const
    { return (bool)_rCell; }

protected:

    LineSegment2D _shape;

    Label _i, _j, _lid;

    const StructuredGrid2D &_grid;

    const Cell *_lCell, *_rCell;
};

#endif
