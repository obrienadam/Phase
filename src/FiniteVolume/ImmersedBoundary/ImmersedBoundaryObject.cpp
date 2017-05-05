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

void ImmersedBoundaryObject::setInternalCells()
{
    updateCells();
    constructStencils();
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
    std::vector<std::pair<Point2D, Scalar>> bps;

    for(const GhostCellStencil& st: stencils_)
        bps.push_back(std::make_pair(st.boundaryPoint(), p(st.cell())));

    std::sort(bps.begin(), bps.end(), [this](const std::pair<Point2D, Scalar>& pt1, const std::pair<Point2D, Scalar>& pt2){
        Vector2D rVec1 = pt1.first - this->shape().centroid();
        Vector2D rVec2 = pt2.first - this->shape().centroid();

        Scalar theta1 = atan2(rVec1.y, rVec1.x);
        Scalar theta2 = atan2(rVec2.y, rVec2.x);

        return (theta1 < 0. ? theta1 + 2.*M_PI: theta1) < (theta2 < 0. ? theta2 + 2.*M_PI: theta2);
    });

    normalForce_ = Vector2D(0., 0.);

    for(int i = 0, end = bps.size(); i < end; ++i)
    {
        const Point2D &a = bps[i].first;
        const Point2D &b = bps[(i + 1)%end].first;

        Scalar pa = bps[i].second;
        Scalar pb = bps[(i + 1)%end].second;

        Scalar pc = (pa + pb)/2;
        Vector2D norm = (b - a).normalVec();
        normalForce_ -= pc*norm;
    }
}

void ImmersedBoundaryObject::computeShearForce(const ScalarFiniteVolumeField &mu, const VectorFiniteVolumeField &u)
{
    shearForce_ = Vector2D(0., 0.);
}

void ImmersedBoundaryObject::updateCells()
{
    CellZone &fluidCells = grid_.cellZone("fluid");
    freshlyClearedCells_->clear();

    for (Cell &cell: *cells_)
        if (!isInIb(cell.centroid()))
        {
            freshlyClearedCells_->push_back(cell); //- Freshly cleared cells need a time step correction
            fluidCells.moveToGroup(cell); //- This will remove the cell from the ib cell zone
        }

    grid_.setCellsActive(*freshlyClearedCells_);

    ibCells_->clear();
    solidCells_->clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            for (Cell &cell: fluidCells.cellCentersWithin(
                    *(Circle *) shapePtr_.get())) //- The circle method is much more efficient
                cells_->moveToGroup(cell);
            break;
        case Shape2D::BOX:
            for (Cell &cell: fluidCells.cellCentersWithin(
                    *(Box *) shapePtr_.get())) //- The box method is much more efficient
                cells_->moveToGroup(cell);
            break;
        default:
            for (Cell &cell: fluidCells.cellCentersWithin(*shapePtr_))
                cells_->moveToGroup(cell);
    }

    auto isIbCell = [this](const Cell &cell) {
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                return true;

        //for (const DiagonalCellLink &dg: cell.diagonals())
        //    if (!isInIb(dg.cell().centroid()))
        //        return true;

        return false;
    };

    for (Cell &cell: *cells_)
        if (isIbCell(cell))
            ibCells_->push_back(cell);
        else
            solidCells_->push_back(cell);

    grid_.setCellsActive(*ibCells_);
    grid_.setCellsInactive(*solidCells_);
    constructStencils();
}

//- Protected methods

void ImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: *ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *shapePtr_, grid_));
}
