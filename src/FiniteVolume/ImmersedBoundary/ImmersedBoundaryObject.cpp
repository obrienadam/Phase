#include <memory>

#include "ImmersedBoundaryObject.h"
#include "LineSegment2D.h"
#include "QuadraticNormalInterpolation.h"

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
    shapePtr_ = std::shared_ptr<Shape2D>(new Circle(center, radius));
}

void ImmersedBoundaryObject::initPolygon(const std::vector<Point2D> &vertices)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Polygon(vertices));
}

std::pair<Point2D, Vector2D> ImmersedBoundaryObject::intersectionStencil(const Point2D &ptA, const Point2D &ptB) const
{
    const Point2D xc = shape().intersections(LineSegment2D(ptA, ptB))[0];
    const std::pair<Point2D, Point2D> edge = shapePtr_->nearestEdge(xc);

    return std::make_pair(
            xc, -(edge.second - edge.first).normalVec()
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
        if (!shapePtr_->isBoundedBy(cell.centroid(), 1e-10))
        {
            freshlyClearedCells_->push_back(cell); //- Freshly cleared cells need a time step correction
            fluidCells.moveToGroup(cell); //- This will remove the cell from the ib cell zone
        }

    grid_.setCellsActive(*freshlyClearedCells_);

    ibCells_->clear();
    solidCells_->clear();

    for (Cell &cell: fluidCells.cellCentersWithin(*shapePtr_, 1e-10))
        cells_->moveToGroup(cell);

    auto isIbCell = [this](const Cell &cell) {
        for (const InteriorLink &nb: cell.neighbours())
            if (!shapePtr_->isBoundedBy(nb.cell().centroid(), 1e-10))
                return true;

        for (const DiagonalCellLink &dg: cell.diagonals())
            if (!shapePtr_->isBoundedBy(dg.cell().centroid(), 1e-10))
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
