#include <memory>

#include "ImmersedBoundaryObject.h"
#include "TranslatingMotion.h"
#include "OscillatingMotion.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string &name,
                                               Label id,
                                               FiniteVolumeGrid2D &grid)
        :
        name_(name),
        grid_(grid),
        id_(id)
{
    cells_ = CellZone("Cells", grid.cellZoneRegistry());

    zoneRegistry_ = std::make_shared<CellZone::ZoneRegistry>();
    ibCells_ = CellZone("IbCells", zoneRegistry_);
    solidCells_ = CellZone("SolidCells", zoneRegistry_);
    freshCells_ = CellZone("FreshCells", zoneRegistry_);
    deadCells_ = CellZone("DeadCells", zoneRegistry_);
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

void ImmersedBoundaryObject::setMotion(std::shared_ptr<Motion> motion)
{
    motion_ = motion;
}

void ImmersedBoundaryObject::setZone(CellZone &zone)
{
    fluid_ = &zone;
    cells_ = CellZone("Cells", zone.registry());
}

void ImmersedBoundaryObject::clear()
{
    fluid_->add(cells_);
    ibCells_.clear();
    solidCells_.clear();
    freshCells_.clear();
    deadCells_.clear();
    grid_.setCellsActive(fluid_->begin(), fluid_->end());
}

LineSegment2D ImmersedBoundaryObject::intersectionLine(const LineSegment2D& ln) const
{
    auto xc = shapePtr_->intersections(ln);
    return xc.empty() ? ln: LineSegment2D(ln.ptA(), xc[0]);
}

Vector2D ImmersedBoundaryObject::nearestEdgeNormal(const Point2D &pt) const
{
    switch(shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            return (shapePtr_->centroid() - pt).unitVec();
        default:
            throw Exception("ImmersedBoundaryObject", "nearestEdgeNormal", "not implemented for specified shape.");
    }
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

Vector2D ImmersedBoundaryObject::acceleration() const
{
    return motion_ ? motion_->acceleration(): Vector2D(0., 0.);
}

Vector2D ImmersedBoundaryObject::acceleration(const Point2D &point) const
{
    return motion_ ? motion_->acceleration(point): Vector2D(0., 0.);
}

Vector2D ImmersedBoundaryObject::velocity() const
{
    return motion_ ? motion_->velocity(): Vector2D(0., 0.);
}

Vector2D ImmersedBoundaryObject::velocity(const Point2D &point) const
{
    return motion_ ? motion_->velocity(point): Vector2D(0., 0.);
}

void ImmersedBoundaryObject::update(Scalar timeStep)
{
    if(motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

void ImmersedBoundaryObject::updateCells()
{
    clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
        {
            auto cells = fluid_->itemsWithin(*std::static_pointer_cast<Circle>(shapePtr_));
            cells_.add(cells.begin(), cells.end());
        }
            break;
        case Shape2D::BOX:
        {
            auto cells = fluid_->itemsWithin(*std::static_pointer_cast<Box>(shapePtr_));
            cells_.add(cells.begin(), cells.end());
        }
            break;
        default:
        {
            auto cells = fluid_->itemsWithin(*shapePtr_);
            cells_.add(cells.begin(), cells.end());
        }
            break;
    }

    solidCells_.add(cells_); //- By default solid cells are not made inactive
}

Equation<Vector2D> ImmersedBoundaryObject::solidVelocity(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    for(const Cell& cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    return eqn;
}

void ImmersedBoundaryObject::clearFreshCells()
{
    fluid_->add(freshCells_); //- Should also clear cells from cells_
    freshCells_.clear();
}
