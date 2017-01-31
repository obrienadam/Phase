#include "ImmersedBoundaryObject.h"
#include "LineSegment2D.h"

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
    const auto &bi = imagePointInterpolation(cell);

    std::vector<Scalar> vals;
    for(const Cell &kCell: ipCells)
        vals.push_back(field[kCell.id()]);

    return bi(vals, ip);
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

void ImmersedBoundaryObject::setInternalCells()
{
    flagIbCells();
    constructStencils();
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
    auto internalCells = grid_.localActiveCells().rangeSearch(shape(), 1e-8);

    std::vector<Label> cells, ibCells, solidCells;

    for(const Cell& cell: internalCells)
        cells.push_back(cell.id());

    grid_.setCellsInactive(cells);

    for(const Cell &cell: internalCells)
    {
        bool isCellActive = false;

        for(const InteriorLink &nb: cell.neighbours())
        {
            if(nb.cell().isActive())
            {
                ibCells.push_back(cell.id());
                isCellActive = true;
                break;
            }
        }

        if(isCellActive)
            continue;

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            if(dg.cell().isActive())
            {
                ibCells.push_back(cell.id());
                isCellActive = true;
                break;
            }
        }

        if(isCellActive)
            continue;

        solidCells.push_back(cell.id());
    }

    grid_.setCellsActive(ibCells);
    cells_ = &grid_.createCellZone(name_ + "IbCells", ibCells);
    grid_.createCellZone(name_ + "SolidCells", solidCells);
}

void ImmersedBoundaryObject::constructStencils()
{
    stencilPoints_.clear();
    imagePointStencils_.clear();

    for(const Cell& cell: cells())
    {
        const Point2D bp = shape().nearestIntersect(cell.centroid()), ip = 2.*bp - cell.centroid();

        stencilPoints_[cell.id()] = std::make_pair(bp, ip);

        const auto& kNN = grid_.findNearestNode(ip).cells();

        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        imagePointStencils_[cell.id()] = std::make_pair(kNN, BilinearInterpolation(centroids));
    }
}
