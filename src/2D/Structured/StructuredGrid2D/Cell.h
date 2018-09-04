#ifndef PHASE_CELL_H
#define PHASE_CELL_H

#include "Orientation.h"
#include "Geometry/Box.h"
#include "System/StaticVector.h"

class StructuredGrid2D;
class Face;

class Cell
{
public:

    Cell(const StructuredGrid2D &grid, Label i, Label j);

    const Box& shape() const
    { return _shape; }

    Point2D centroid() const
    { return _shape.centroid(); }

    Label i() const
    { return _i; }

    Label j() const
    { return _j; }

    Label id() const
    { return _lid; }

    Label lid() const
    { return _lid; }

    Label gid() const
    { return _gid; }

    const StructuredGrid2D &grid() const
    { return _grid; }

    const Cell &cell(Orientation dir, int offset) const;

    Scalar dh(Orientation dir, int offset) const;

    const Face &face(Orientation dir) const;

protected:

    Box _shape;

    Label _i, _j, _lid, _gid;    //- Stencils

    const StructuredGrid2D &_grid;

};

#endif
