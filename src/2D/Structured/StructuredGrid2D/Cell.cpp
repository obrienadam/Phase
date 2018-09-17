#include "Cell.h"
#include "StructuredGrid2D.h"

Cell::Cell(const StructuredGrid2D &grid, Label i, Label j)
    :
      _grid(grid),
      _i(i),
      _j(j),
      _lid(grid.cells().size()),
      _gid(grid.cells().size())
{
    _shape = Box(grid.node(_i, _j), grid.node(_i + 1, _j + 1));
}

Cell::Cell(const StructuredGrid2D &grid, Label i, Label j, Label gid)
    :
      Cell(grid, i, j)
{
    _gid = gid;
}

void Cell::initStencils()
{
    for(auto zeta: Coordinates::DIRECTIONS)
        if(_grid.maxInc(*this, zeta) > 0)
            _interiorFaceStencils.push_back(InteriorFaceStencil(*this, zeta));
}
