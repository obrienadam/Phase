#include <memory>

#include "ImmersedBoundaryObject.h"
#include "LineSegment2D.h"
#include "QuadraticNormalInterpolation.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               Label id,
                                               FiniteVolumeGrid2D &grid)
    :
      name_(name),
      grid_(grid),
      id_(id)
{
    freshlyClearedCells_ = &grid.createCellGroup(name_ + "FreshlyClearedCells", std::vector<Label>());
}

void ImmersedBoundaryObject::initCircle(const Point2D &center, Scalar radius)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Circle(center, radius));
}

void ImmersedBoundaryObject::initPolygon(const std::vector<Point2D> &vertices)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Polygon(vertices));
}

std::pair<Point2D, Vector2D> ImmersedBoundaryObject::intersectionStencil(const Point2D& ptA, const Point2D& ptB) const
{
    const LineSegment2D line(ptA, ptB);
    const Point2D xc = intersection(line, shape())[0];

    const std::pair<Point2D, Point2D> edge = shapePtr_->nearestEdge(xc);

    return std::make_pair(
                xc, -(edge.second - edge.first).normalVec()
                );
}

void ImmersedBoundaryObject::setInternalCells()
{
    flagIbCells();
    constructStencils();
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

template<>
Scalar ImmersedBoundaryObject::getBoundaryRefValue<Scalar>(const std::string& name) const
{
    return boundaryRefScalars_.find(name)->second;
}

template<>
Vector2D ImmersedBoundaryObject::getBoundaryRefValue<Vector2D>(const std::string& name) const
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
    catch(...)
    {
        boundaryRefVectors_[name] = Vector2D(value);
    }
}

void ImmersedBoundaryObject::updateCells()
{
    ibCells_->clear();
    solidCells_->clear();

    CellZone& fluidCells = grid_.cellZone("fluid");

    for(Cell& cell: cells_->cells())
    {
        fluidCells.moveToGroup(cell);
        grid_.setCellActive(cell);
    }

    flagIbCells();
    constructStencils();
}

//- Protected methods

void ImmersedBoundaryObject::flagIbCells()
{
    auto internalCells = grid_.localActiveCells().rangeSearch(shape(), 1e-10);

    std::vector<Label> cells, ibCells, solidCells;

    for(const Cell& cell: internalCells)
        cells.push_back(cell.id());

    //- Create a new cell zone for the IB (this does not affect if the cell is active or not)
    cells_ = &grid_.createCellZone(name_ + "Cells", cells);

    //- Sort the new zone into ib cells and solid cells
    for(const Cell &cell: *cells_)
    {
        bool isIbCell = false;

        for(const InteriorLink &nb: cell.neighbours())
        {
            if(!shape().isCovered(nb.cell().centroid()))
            {
                ibCells.push_back(cell.id());
                isIbCell = true;
                break;
            }
        }

        if(isIbCell)
            continue;

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            if(!shape().isCovered(dg.cell().centroid()))
            {
                ibCells.push_back(cell.id());
                isIbCell = true;
                break;
            }
        }

        if(isIbCell)
            continue;

        solidCells.push_back(cell.id());
    }

    grid_.setCellsInactive(solidCells);
    ibCells_ = &grid_.createCellGroup(name_ + "IbCells", ibCells);
    solidCells_ = &grid_.createCellGroup(name_ + "SolidCells", solidCells);
}

void ImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for(const Cell& cell: *ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *shapePtr_, grid_));
}
