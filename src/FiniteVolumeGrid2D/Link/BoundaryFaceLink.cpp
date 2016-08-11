#include "BoundaryFaceLink.h"
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

BoundaryLink::BoundaryLink(const BoundaryLink &other)
    :
      BoundaryLink(other.self_, other.face_)
{

}

const Face& BoundaryLink::face() const
{
    return face_;
}
