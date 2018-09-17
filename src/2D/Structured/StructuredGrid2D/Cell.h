#ifndef PHASE_CELL_H
#define PHASE_CELL_H

#include "Coordinates.h"
#include "Geometry/Box.h"
#include "System/StaticVector.h"

#include "InteriorFaceStencil.h"

class StructuredGrid2D;
class Face;

class Cell
{
public:

    Cell(const StructuredGrid2D &grid, Label i, Label j);

    Cell(const StructuredGrid2D &grid, Label i, Label j, Label gid);

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

    void setgid(Label gid)
    { _gid = gid; }

    const StructuredGrid2D &grid() const
    { return _grid; }

    //- Stencils
    void initStencils();

    const std::vector<InteriorFaceStencil> &interiorFaceStencils() const
    { return _interiorFaceStencils; }

protected:

    std::vector<InteriorFaceStencil> _interiorFaceStencils;

    Box _shape;

    Label _i, _j, _lid, _gid;

    const StructuredGrid2D &_grid;

};

#endif
