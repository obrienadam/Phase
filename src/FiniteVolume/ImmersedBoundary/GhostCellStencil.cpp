#include "GhostCellStencil.h"
#include "FiniteVolumeGrid2D.h"

GhostCellStencil::GhostCellStencil(const Cell &cell, const Shape2D &shape, const FiniteVolumeGrid2D &grid)
        :
          ImmersedBoundaryStencil(cell)
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

Vector2D GhostCellStencil::ipGrad(const ScalarFiniteVolumeField &field) const
{
    Matrix mat(4, 3), b(4, 1);

    int i = 0;
    Scalar ipVal = ipValue(field);
    for (const Cell &cell: ipCells_)
    {
        mat(i, 0) = cell.centroid().x - ip_.x;
        mat(i, 1) = cell.centroid().y - ip_.y;
        b(i, 0) = field(cell) - ipVal;
        ++i;
    }

    mat.solve(b);
    return Vector2D(b(0, 0), b(0, 1));
}
