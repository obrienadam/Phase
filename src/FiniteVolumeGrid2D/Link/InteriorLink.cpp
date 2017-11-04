#include "InteriorLink.h"
#include "Cell.h"

//- Interior link

InteriorLink::InteriorLink(const Cell &self, const Face &face, const Cell &cell)
        :
        CellLink(self, cell),
        face_(face)
{
    outwardNorm_ = face_.outwardNorm(self_.centroid());
    rFaceVec_ = face_.centroid() - self_.centroid();
}

InteriorLink::InteriorLink(const InteriorLink &other)
        :
        InteriorLink(other.self_, other.face_, other.cell_)
{

}

Scalar InteriorLink::volumeWeight() const
{
    return cell_.volume() / (self_.volume() + cell_.volume());
}


Scalar InteriorLink::distanceWeight() const
{
    Scalar l1 = (face_.centroid() - self_.centroid()).mag();
    Scalar l2 = (face_.centroid() - cell_.centroid()).mag();
    return l2 / (l1 + l2);
}

Scalar InteriorLink::distanceSqrWeight() const
{
    Scalar l1 = (face_.centroid() - self_.centroid()).magSqr();
    Scalar l2 = (face_.centroid() - cell_.centroid()).magSqr();
    return l2 / (l1 + l2);
}
