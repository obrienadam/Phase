#include "InteriorFaceLink.h"
#include "Cell.h"

//- Interior link

InteriorLink::InteriorLink(const Cell &self, const Face &face, const Cell &cell)
    :
      BoundaryLink(self, face),
      cell_(cell)
{
    rCellVec_ = cell.centroid() - self.centroid();
}
