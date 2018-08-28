#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "Cell.h"

class Face
{
public:

    enum Orientation{I, J};

    Face(const StructuredGrid2D &grid, Orientation orientation, Label i, Label j);

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
