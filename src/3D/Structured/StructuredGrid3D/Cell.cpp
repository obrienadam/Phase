#include "Cell.h"
#include "StructuredGrid3D.h"

Cell::Cell(const StructuredGrid3D &grid, Label i, Label j, Label k)
    :
      _grid(grid),
      _i(i),
      _j(j),
      _k(k),
      _id(grid.cells().size())
{
    _shape = RectangularPrism(grid.node(i, j, k), grid.node(i + 1, j + 1, k + 1));
}

const Cell &Cell::nb(Cell::Index idx, int offset) const
{
    switch(idx)
    {
    case I:
        return _grid(_i + offset, _j, _k);
    case J:
        return _grid(_i, _j + offset, _k);
    case K:
        return _grid(_i, _j, _k + offset);
    }
}

void Cell::initStencils(int order)
{

}
