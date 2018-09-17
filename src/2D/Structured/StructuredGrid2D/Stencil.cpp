#include "Stencil.h"
#include "Cell.h"
#include "StructuredGrid2D.h"

Stencil::Stencil(const Cell &cell, Coordinates::Direction zeta)
    :
      _cell(cell),
      _zeta(zeta)
{

}

const Cell &Stencil::cell(int inc) const
{
    return _cell.grid().cell(_cell, _zeta, inc);
}

Vector2D Stencil::rc(int inc) const
{
    return cell(inc).centroid() - _cell.centroid();
}
