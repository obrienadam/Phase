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

InteriorLink::InteriorLink(const InteriorLink &other)
        :
        InteriorLink(other.self_, other.face_, other.cell_)
{

}

void InteriorLink::init()
{
    BoundaryLink::init();
    rCellVec_ = cell_.centroid() - self_.centroid();
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

const Cell &InteriorLink::cell() const
{
    return cell_;
}
