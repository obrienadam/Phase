#include "BoundaryLink.h"
#include "Face.h"
#include "Cell.h"

BoundaryLink::BoundaryLink(const Cell &self, const Face &face)
    :
      Link(self),
      face_(face)
{
    rFaceVec_ = face_.centroid() - self_.centroid();
    outwardNorm_ = face_.outwardNorm(self_.centroid());
}

BoundaryLink::BoundaryLink(const BoundaryLink &other)
    :
      BoundaryLink(other.self_, other.face_)
{

}

const Face &BoundaryLink::face() const
{
    return face_;
}
