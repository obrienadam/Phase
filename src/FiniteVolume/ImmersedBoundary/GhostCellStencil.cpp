#include "GhostCellStencil.h"
#include "FiniteVolumeGrid2D.h"

GhostCellStencil::GhostCellStencil(const Cell &cell, const Shape2D &shape, const FiniteVolumeGrid2D &grid)
        :
        cell_(cell)
{
    bp_ = shape.nearestIntersect(cell_.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    ipCells_ = grid.findNearestNode(ip_).cells();

    std::vector<Point2D> centroids;
    for (const Cell &cell: ipCells_)
        centroids.push_back(cell.centroid());

    interpolator_ = std::unique_ptr<Interpolation>(new BilinearInterpolation(centroids));

    ipCoeffs_ = (*interpolator_)(ip_);
}

Scalar GhostCellStencil::ipValue(const ScalarFiniteVolumeField &field) const
{
    std::vector<Scalar> vals;

    for (const Cell &cell: ipCells_)
        vals.push_back(field(cell));

    return (*interpolator_)(vals, ip_);
}

Vector2D GhostCellStencil::ipValue(const VectorFiniteVolumeField &field) const
{
    std::vector<Vector2D> vals;

    for (const Cell &cell: ipCells_)
        vals.push_back(field(cell));

    return (*interpolator_)(vals, ip_);
}
