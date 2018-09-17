#include "InteriorFaceStencil.h"
#include "Cell.h"
#include "Face.h"
#include "StructuredGrid2D.h"

InteriorFaceStencil::InteriorFaceStencil(const Cell &cell, Coordinates::Direction zeta)
    :
      Stencil(cell, zeta),
      _face(cell.grid().face(cell, zeta))
{
    _rf = _face.centroid() - _cell.centroid();
    _sf = _face.sf(cell.centroid());
}
