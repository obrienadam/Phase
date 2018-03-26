#include "ImmersedBoundaryObject.h"
#include "ImmersedBoundary.h"
#include "NotImplementedException.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string &name,
                                               Label id,
                                               const ImmersedBoundary &ib,
                                               const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        name_(name),
        ib_(&ib),
        grid_(grid),
        id_(id)
{
    cells_ = CellZone("Cells", grid->cellZoneRegistry());
    zoneRegistry_ = std::make_shared<CellZone::ZoneRegistry>();
    ibCells_ = CellZone("IbCells", zoneRegistry_);
    solidCells_ = CellZone("SolidCells", zoneRegistry_);
    freshCells_ = CellZone("FreshCells", zoneRegistry_);

    force_ = Vector2D(0., 0.);
    torque_ = 0.;
}

void ImmersedBoundaryObject::initCircle(const Point2D &center, Scalar radius)
{
    shape_ = std::unique_ptr<Circle>(new Circle(center, radius));
}


void ImmersedBoundaryObject::initBox(const Point2D &center, Scalar width, Scalar height)
{
    shape_ = std::unique_ptr<Box>(new Box(
            Point2D(center.x - width / 2., center.y - height / 2.),
            Point2D(center.x + width / 2., center.y + height / 2.)
    ));
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
}

LineSegment2D ImmersedBoundaryObject::intersectionLine(const LineSegment2D &ln) const
{
    auto xc = shape_->intersections(ln);

    if (xc.empty())
    {
        Point2D pts[] = {ln.ptA(), ln.ptB()};
        xc.push_back(shape_->closest(pts, pts + 2));
    }

    return LineSegment2D(ln.ptA(), xc[0]);
}

LineSegment2D ImmersedBoundaryObject::intersectionLine(const Point2D &ptA, const Point2D &ptB) const
{
    return intersectionLine(LineSegment2D(ptA, ptB));
}

Vector2D ImmersedBoundaryObject::nearestEdgeNormal(const Point2D &pt) const
{
    switch (shape_->type())
    {
        case Shape2D::CIRCLE:return (shape_->centroid() - pt).unitVec();
        case Shape2D::BOX:
        case Shape2D::POLYGON:
        {
            auto edge = shape_->nearestEdge(pt);
            return dot(edge.norm(), shape_->centroid() - edge.center()) > 0. ? edge.norm().unitVec()
                                                                             : -edge.norm().unitVec();
        }
        default:throw Exception("ImmersedBoundaryObject", "nearestEdgeNormal", "not implemented for specified shape.");
    }
}

void ImmersedBoundaryObject::addBoundary(const std::string &name, BoundaryType bType, Scalar ref)
{
    boundaryTypes_[name] = bType;
    boundaryRefScalars_[name] = ref;
}

void ImmersedBoundaryObject::addBoundary(const std::string &name, BoundaryType bType, const Vector2D &ref)
{
    boundaryTypes_[name] = bType;
    boundaryRefVectors_[name] = ref;
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, const std::string &boundaryType)
{
    if (boundaryType == "fixed")
        addBoundaryType(name, FIXED);
    else if (boundaryType == "normal_gradient")
        addBoundaryType(name, NORMAL_GRADIENT);
    else if (boundaryType == "partial_slip")
        addBoundaryType(name, PARTIAL_SLIP);
    else
        throw Exception("ImmersedBoundaryObject",
                        "ImmersedBoundaryObject",
                        "invalid boundary type \"" + boundaryType + "\".");
}

template<>
Scalar ImmersedBoundaryObject::getBoundaryRefValue<Scalar>(const std::string &name) const
{
    auto it = boundaryRefScalars_.find(name);
    return it == boundaryRefScalars_.end() ? 0. : it->second;
}

template<>
Vector2D ImmersedBoundaryObject::getBoundaryRefValue<Vector2D>(const std::string &name) const
{
    auto it = boundaryRefVectors_.find(name);
    return it == boundaryRefVectors_.end() ? Vector2D(0., 0.) : it->second;
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

Vector2D ImmersedBoundaryObject::velocity(const Point2D &point) const
{
    return motion_ ? motion_->velocity(point) : Vector2D(0., 0.);
}

Vector2D ImmersedBoundaryObject::acceleration(const Point2D &point) const
{
    return motion_ ? motion_->acceleration(point) : Vector2D(0., 0.);
}

Scalar ImmersedBoundaryObject::theta() const
{
    return motion_ ? motion_->theta() : 0.;
}

Scalar ImmersedBoundaryObject::omega() const
{
    return motion_ ? motion_->omega() : 0.;
}

Scalar ImmersedBoundaryObject::alpha() const
{
    return motion_ ? motion_->alpha() : 0.;
}

std::vector<Point2D> ImmersedBoundaryObject::forcingPoints() const
{
    std::vector<Point2D> points(ibCells_.size());

    std::transform(ibCells_.begin(), ibCells_.end(), points.begin(), [this](const Cell &cell)
    {
        return nearestIntersect(cell.centroid());
    });

    std::sort(points.begin(), points.end(), [this](const Point2D &lhs, const Point2D &rhs)
    {
        return (lhs - shape_->centroid()).angle() < (rhs - shape_->centroid()).angle();
    });

    return points;
}

void ImmersedBoundaryObject::computeForce(Scalar rho,
                                          Scalar mu,
                                          const VectorFiniteVolumeField &u,
                                          const ScalarFiniteVolumeField &p,
                                          const Vector2D &g)
{
    throw NotImplementedException("ImmersedBoundaryObject", "computeForce");
}

void ImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                          const ScalarFiniteVolumeField &mu,
                                          const VectorFiniteVolumeField &u,
                                          const ScalarFiniteVolumeField &p,
                                          const Vector2D &g)
{
    throw NotImplementedException("ImmersedBoundaryObject", "computeForce");
}

void ImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                          const ScalarFiniteVolumeField &mu,
                                          const VectorFiniteVolumeField &u,
                                          const ScalarFiniteVolumeField &p,
                                          const ScalarFiniteVolumeField &gamma,
                                          const SurfaceTensionForce &ft,
                                          const Vector2D &g)
{
    throw NotImplementedException("ImmersedBoundaryObject", "computeForce");
}

void ImmersedBoundaryObject::update(Scalar timeStep)
{
    if (motion_)
    {
        motion_->update(timeStep);
        shape_->move(motion_->position());
        updateCells();
    }
}

void ImmersedBoundaryObject::updateCells()
{
    clear();
    auto cells = fluid_->itemsWithin(*shape_);
    cells_.add(cells.begin(), cells.end());
    solidCells_.add(cells_); //- By default solid cells are not made inactive
}

Equation<Vector2D> ImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    return eqn;
}

Equation<Scalar> ImmersedBoundaryObject::pressureBcs(ScalarFiniteVolumeField &p) const
{
    throw NotImplementedException("ImmersedBoundaryObject", "pressureBcs");
}

void ImmersedBoundaryObject::computeBoundaryForcing(const VectorFiniteVolumeField& u,
                                                    Scalar timeStep,
                                                    VectorFiniteVolumeField &fb) const
{
    throw NotImplementedException("ImmersedBoundaryObject", "computeBoundaryForcing");
}

void ImmersedBoundaryObject::clearFreshCells()
{
    fluid_->add(freshCells_); //- Should also clear cells from cells_
    freshCells_.clear();
}
