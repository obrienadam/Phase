#include "Link.h"
#include "Cell.h"
#include "Face.h"

//- Boundary link

BoundaryLink::BoundaryLink(const Cell& self, const Face &face)
    :
      Link(self),
      face_(face)
{
    rFaceVec_ = face.centroid() - self.centroid();
    outwardNorm_ = face.outwardNorm(self.centroid());
}

//- Interior link

InteriorLink::InteriorLink(const Cell &self, const Face &face, const Cell &cell, Direction direction)
    :
      BoundaryLink(self, face),
      cell_(cell),
      direction_(direction)
{
    rCellVec_ = cell.centroid() - self.centroid();
}
