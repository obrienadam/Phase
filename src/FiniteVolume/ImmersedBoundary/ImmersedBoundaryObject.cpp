#include "ImmersedBoundaryObject.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const FiniteVolumeGrid2D& grid,
                                               const Point2D& center,
                                               Scalar radius)
    :
      name_(name),
      grid_(grid)
{
    shapePtr_ = std::shared_ptr<Shape2D>(new Circle(center, radius));
}

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const FiniteVolumeGrid2D& grid,
                                               const std::vector<Point2D> &vertices)
    :
      name_(name),
      grid_(grid)
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
    const Line2D line(ptA, (ptB - ptA).normalVec());
    const Point2D xc = nearestPoint(ptA, shapePtr_->intersections(line));
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
    auto internalCells = grid_.activeCells().rangeSearch(shape());
    std::vector<Label> internalCellIds;

    for(const Cell &cell: internalCells)
        internalCellIds.push_back(cell.id());

    grid_.moveCellsToInactiveCellGroup(internalCellIds);
    internalCellIds.clear();

    for(const Cell &cell: internalCells)
    {
        bool isCellActive = false;

        for(const InteriorLink &nb: cell.neighbours())
        {
            if(nb.cell().isActive())
            {
                internalCellIds.push_back(cell.id());
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
                internalCellIds.push_back(cell.id());
                break;
            }
        }
    }

    grid_.moveCellsToCellGroup(name_ + "_cells", internalCellIds);
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
