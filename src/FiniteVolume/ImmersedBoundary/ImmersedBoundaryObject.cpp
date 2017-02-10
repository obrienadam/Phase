#include <memory>

#include "ImmersedBoundaryObject.h"
#include "LineSegment2D.h"
#include "QuadraticNormalInterpolation.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const Point2D& center,
                                               Scalar radius,
                                               Label id,
                                               FiniteVolumeGrid2D &grid)
    :
      name_(name),
      grid_(grid),
      id_(id)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Circle(center, radius));
}

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const std::vector<Point2D> &vertices,
                                               Label id,
                                               FiniteVolumeGrid2D& grid)
    :
      name_(name),
      grid_(grid),
      id_(id)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Polygon(vertices));
}

Scalar ImmersedBoundaryObject::imagePointVal(const Cell &cell, const ScalarFiniteVolumeField& field) const
{
    const Point2D &ip = imagePoint(cell);
    const auto &ipCells = imagePointCells(cell);
    const auto &interpolator = imagePointInterpolation(cell);

    std::vector<Scalar> vals;
    for(const Cell &kCell: ipCells)
        vals.push_back(field[kCell.id()]);

    return interpolator(vals, ip);
}

Vector2D ImmersedBoundaryObject::imagePointVal(const Cell &cell, const VectorFiniteVolumeField& field) const
{
    const Point2D &ip = imagePoint(cell);
    const auto &ipCells = imagePointCells(cell);
    const auto &bi = imagePointInterpolation(cell);

    std::vector<Vector2D> vals;
    for(const Cell &kCell: ipCells)
        vals.push_back(field[kCell.id()]);

    return bi(vals, ip);
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

void ImmersedBoundaryObject::setInternalCells(const Communicator &comm)
{
    flagIbCells();
    constructStencils();
    comm.printf("IB cells set for object \"%s\".\n", name().c_str());
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

void ImmersedBoundaryObject::addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue)
{
    boundaryRefValues_[name] = boundaryRefValue;
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
    stencilPoints_.clear();
    imagePointStencils_.clear();

    for(const Cell& cell: ibCells())
    {
        const Point2D bp = shape().nearestIntersect(cell.centroid()), ip = 2.*bp - cell.centroid();

        stencilPoints_[cell.id()] = std::make_pair(bp, ip);

        const auto& kNN = grid_.findNearestNode(ip).cells();

        for(const Cell& cell: kNN)
            if(solidCells_->isInGroup(cell))
                throw Exception("ImmersedBoundaryActive", "constructStencils", "attempted to add solid cell to stencil.");

        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        switch(interpolationType_)
        {
        case BILINEAR:
            imagePointStencils_[cell.id()] = std::make_pair(kNN, std::unique_ptr<Interpolation>(new BilinearInterpolation(centroids)));
            break;
        case QUADRATIC_NORMAL:
            imagePointStencils_[cell.id()] = std::make_pair(kNN, std::unique_ptr<Interpolation>(new QuadraticNormalInterpolation(centroids, ip - bp)));
            break;
        }
    }
}
