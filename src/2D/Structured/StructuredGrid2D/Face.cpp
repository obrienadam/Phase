#include "Face.h"
#include "StructuredGrid2D.h"

Face::Face(const StructuredGrid2D &grid, Coordinate coord, Label i, Label j)
    :
      _i(i),
      _j(j),
      _lid(grid.ifaces().size() + grid.jfaces().size()),
      _grid(grid)
{
    switch (coord)
    {
    case I:
        _shape = LineSegment2D(_grid.node(_i, _j), _grid.node(_i, _j + 1));
        break;
    case J:
        _shape = LineSegment2D(_grid.node(_i, _j), _grid.node(_i + 1, _j));
    }
}
