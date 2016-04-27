#include "BoundaryLink.h"
#include "Face.h"
#include "Cell.h"

BoundaryLink::BoundaryLink(const Cell& self, const Face &face)
    :
      Link(self),
      face_(face)
{
    rFaceVec_ = face.centroid() - self.centroid();
    outwardNorm_ = face.outwardNorm(self.centroid());
}
