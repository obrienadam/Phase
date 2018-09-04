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

const Cell &Cell::cell(Orientation dir, int offset) const
{
    return _grid.cell(*this, dir, offset);
}

Scalar Cell::dh(Orientation dir, int offset) const
{
    return _grid.dh(*this, dir, offset);
}

const Face &Cell::face(Orientation dir) const
{
    return _grid.face(*this, dir);
}
