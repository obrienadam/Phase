#include <memory>

#include "ImmersedBoundaryObject.h"
#include "Box.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string &name,
                                               Label id,
                                               FiniteVolumeGrid2D &grid)
    :
      name_(name),
      grid_(grid),
      id_(id)
{
    cells_ = &grid_.createCellZone(name_ + "Cells");
    ibCells_ = &grid_.createCellGroup(name_ + "IbCells");
    solidCells_ = &grid_.createCellGroup(name_ + "SolidCells");
    freshlyClearedCells_ = &grid.createCellGroup(name_ + "FreshlyClearedCells");
}

void ImmersedBoundaryObject::initCircle(const Point2D &center, Scalar radius)
{
    shapePtr_ = std::shared_ptr<Circle>(new Circle(center, radius));
}


void ImmersedBoundaryObject::initBox(const Point2D &center, Scalar width, Scalar height)
{
    shapePtr_ = std::shared_ptr<Box>(new Box(
                                         Point2D(center.x - width / 2., center.y - height / 2.),
                                         Point2D(center.x + width / 2., center.y + height / 2.)
                                         ));
}

std::pair<Point2D, Vector2D> ImmersedBoundaryObject::intersectionStencil(const Point2D &ptA, const Point2D &ptB) const
{
    auto intersections = shape().intersections(LineSegment2D(ptA, ptB));

    Point2D xc;
    if (intersections.empty()) //- fail safe, in case a point is on an ib
    {
        Point2D nPtA = shape().nearestIntersect(ptA);
        Point2D nPtB = shape().nearestIntersect(ptB);

        if ((nPtA - ptA).magSqr() < (nPtB - ptB).magSqr())
            xc = ptA;
        else
            xc = ptB;
    }
    else
        xc = intersections[0];

    LineSegment2D edge = shapePtr_->nearestEdge(xc);

    return std::make_pair(
                xc, -(edge.ptB() - edge.ptA()).normalVec()
                );
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

template<>
Scalar ImmersedBoundaryObject::getBoundaryRefValue<Scalar>(const std::string &name) const
{
    return boundaryRefScalars_.find(name)->second;
}

template<>
Vector2D ImmersedBoundaryObject::getBoundaryRefValue<Vector2D>(const std::string &name) const
{
    return boundaryRefVectors_.find(name)->second;
}

void ImmersedBoundaryObject::addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue)
{
    boundaryRefScalars_[name] = boundaryRefValue;
}

void ImmersedBoundaryObject::addBoundaryRefValue(const std::string &name, const Vector2D &boundaryRefValue)
{
    boundaryRefVectors_[name] = boundaryRefValue;
}

void ImmersedBoundaryObject::addBoundaryRefValue(const std::string &name, const std::string &value)
{
    try
    {
        boundaryRefScalars_[name] = std::stod(value);
        return;
    }
    catch (...)
    {
        boundaryRefVectors_[name] = Vector2D(value);
    }
}

void ImmersedBoundaryObject::computeNormalForce(const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &p)
{

}

void ImmersedBoundaryObject::computeShearForce(const ScalarFiniteVolumeField &mu, const VectorFiniteVolumeField &u)
{

}
