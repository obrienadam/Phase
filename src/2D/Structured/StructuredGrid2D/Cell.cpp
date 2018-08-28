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
