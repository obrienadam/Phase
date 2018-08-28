#ifndef PHASE_CELL_H
#define PHASE_CELL_H

#include "Geometry/Box.h"

class StructuredGrid2D;

class Cell
{
public:

    enum Face{I_POS, I_NEG, J_POS, J_NEG};

    Cell(const StructuredGrid2D &grid, Label i, Label j);

    const Box& shape() const
    { return _shape; }

    Point2D centroid() const
    { return _shape.centroid(); }

    Label lid() const
    { return _lid; }

    const StructuredGrid2D &grid() const
    { return _grid; }

protected:

    Box _shape;

    Label _i, _j, _lid, _gid;

    const StructuredGrid2D &_grid;
};

#endif
