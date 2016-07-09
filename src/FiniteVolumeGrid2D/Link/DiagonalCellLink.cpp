#include "DiagonalCellLink.h"
#include "Cell.h"

DiagonalCellLink::DiagonalCellLink(const Cell &self, const Cell &cell)
    :
      Link(self),
      cell_(cell)
{
    rCellVec_ = cell.centroid() - self.centroid();
}
