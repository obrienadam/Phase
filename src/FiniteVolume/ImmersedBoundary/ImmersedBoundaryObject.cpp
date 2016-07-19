#include "ImmersedBoundaryObject.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const FiniteVolumeGrid2D& grid,
                                               const Point2D& center,
                                               Scalar radius)
    :
      Circle(center, radius),
      name_(name),
      grid_(grid)
{

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

const std::vector< Ref<const Cell> > ImmersedBoundaryObject::boundingCells(const Point2D &pt) const
{
    const Node &node = grid_.findNearestNode(pt);
    return grid_.activeCells().kNearestNeighbourSearch(node, 4);
}

void ImmersedBoundaryObject::flagIbCells()
{
    auto internalCells = grid_.activeCells().rangeSearch(*this);
    std::vector<size_t> internalCellIds;

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
        const Point2D bp = nearestIntersect(cell.centroid()), ip = 2.*bp - cell.centroid();

        stencilPoints_[cell.id()] = std::make_pair(bp, ip);

        auto kNN = boundingCells(ip);

        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        imagePointStencils_[cell.id()] = std::make_pair(kNN, BilinearInterpolation(centroids));
    }
}
