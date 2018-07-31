#include "GhostCellImmersedBoundaryFixedBcStencil.h"

GhostCellImmersedBoundary::FixedBcStencil::FixedBcStencil(const Cell &cell,
                                                          const ImmersedBoundaryObject &ibObj,
                                                          const FiniteVolumeGrid2D &grid)
    :
      BcStencil(cell, ibObj, grid)
{
    initInterpolationCoeffs(ibObj, grid);
}

GhostCellImmersedBoundary::FixedBcStencil::FixedBcStencil(const Cell &cell,
                                                          const ImmersedBoundaryObject &ibObj,
                                                          const FiniteVolumeGrid2D &grid,
                                                          const Vector2D &r)
    :
      BcStencil(cell, ibObj, grid)
{
    initInterpolationCoeffs(ibObj, grid);
}

void GhostCellImmersedBoundary::FixedBcStencil::initInterpolationCoeffs(const ImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid)
{
    cells_.push_back(std::cref(cell_));
    coeffs_ = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Ainv_ / 2.;
    coeffs_.push_back(0.5);
}
