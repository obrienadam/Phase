#include "GhostCellImmersedBoundaryBcStencil.h"

GhostCellImmersedBoundary::BcStencil::BcStencil(const Cell &cell,
                                                const ImmersedBoundaryObject &ibObj,
                                                const FiniteVolumeGrid2D &grid)
    :
      cell_(std::cref(cell))
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = ibObj.isInIb(cell) ? 2. * bp_ - cell.centroid() : 2. * cell.centroid() - bp_;
    nw_ = ibObj.nearestEdgeUnitNormal(bp_);
    initCoeffMatrix(grid);
}

GhostCellImmersedBoundary::BcStencil::BcStencil(const Cell &cell,
                                                const ImmersedBoundaryObject &ibObj,
                                                const FiniteVolumeGrid2D &grid,
                                                const Vector2D &r)
    :
      cell_(std::cref(cell))
{
    bp_ = ibObj.shape().intersections(Ray2D(cell.centroid(), r))[0];
    ip_ = ibObj.isInIb(cell) ? 2. * bp_ - cell.centroid() : 2. * cell.centroid() - bp_;
    nw_ = ibObj.nearestEdgeUnitNormal(bp_);
    initCoeffMatrix(grid);
}

Scalar GhostCellImmersedBoundary::BcStencil::bpValue(const ScalarFiniteVolumeField &phi) const
{
    return (phi(cell_) + ipValue(phi)) / 2.;
}

Scalar GhostCellImmersedBoundary::BcStencil::ipValue(const ScalarFiniteVolumeField &phi) const
{
    auto val = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Ainv_ * StaticMatrix<4, 1>({
                                                                                                      phi(cells_[0]),
                                                                                                      phi(cells_[1]),
                                                                                                      phi(cells_[2]),
                                                                                                      phi(cells_[3])
                                                                                                  });
    return val(0, 0);
}

Tensor2D GhostCellImmersedBoundary::BcStencil::bpGrad(const VectorFiniteVolumeField &u) const
{
    auto cells = u.grid()->findNearestNode(bp_).cells();
    Point2D x[] = {
        cells[0].get().centroid(),
        cells[1].get().centroid(),
        cells[2].get().centroid(),
        cells[3].get().centroid()
    };

    auto A = StaticMatrix<4, 4>({
                                    x[0].x * x[0].y, x[0].x, x[0].y, 1.,
                                    x[1].x * x[1].y, x[1].x, x[1].y, 1.,
                                    x[2].x * x[2].y, x[2].x, x[2].y, 1.,
                                    x[3].x * x[3].y, x[3].x, x[3].y, 1.
                                });

    auto b = StaticMatrix<4, 2>({
                                    u(cells[0]).x, u(cells[0]).y,
                                    u(cells[1]).x, u(cells[1]).y,
                                    u(cells[2]).x, u(cells[2]).y,
                                    u(cells[3]).x, u(cells[3]).y
                                });

    auto grad = StaticMatrix<2, 4>({
                                       bp_.y, 1., 0., 0.,
                                       bp_.x, 0., 1., 0.
                                   }) * solve(A, b);

    return Tensor2D(grad(0, 0), grad(0, 1), grad(1, 0), grad(1, 1));
}

void GhostCellImmersedBoundary::BcStencil::initCoeffMatrix(const FiniteVolumeGrid2D &grid)
{
    cells_ = grid.findNearestNode(ip_).cells();

    Point2D x[] = {
        cells_[0].get().centroid(),
        cells_[1].get().centroid(),
        cells_[2].get().centroid(),
        cells_[3].get().centroid()
    };

    Ainv_ = inverse<4>({
                           x[0].x * x[0].y, x[0].x, x[0].y, 1.,
                           x[1].x * x[1].y, x[1].x, x[1].y, 1.,
                           x[2].x * x[2].y, x[2].x, x[2].y, 1.,
                           x[3].x * x[3].y, x[3].x, x[3].y, 1.
                       });
}
