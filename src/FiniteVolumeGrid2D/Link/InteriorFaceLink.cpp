#include "InteriorFaceLink.h"
#include "Cell.h"

//- Interior link

InteriorLink::InteriorLink(const Cell &self, const Face &face, const Cell &cell, Direction direction)
    :
      BoundaryLink(self, face),
      cell_(cell),
      direction_(direction)
{
    rCellVec_ = cell.centroid() - self.centroid();
}
