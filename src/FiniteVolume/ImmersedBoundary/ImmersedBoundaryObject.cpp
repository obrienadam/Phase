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
    cutCells_ = &grid.createCellGroup(name_ + "CutCells");
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

        for (const DiagonalCellLink &dg: cell.diagonals())
            if (!isInIb(dg.cell().centroid()))
                return true;

        return false;
    };

    for (Cell &cell: *cells_)
        if (isIbCell(cell))
            ibCells_->push_back(cell);
        else
            solidCells_->push_back(cell);

    //- Flag cut cells, used for some IB methods
    cutCells_->clear();
    for (Cell &cell: fluidCells.cellsOverlapping(shape()))
        cutCells_->push_back(cell);

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
