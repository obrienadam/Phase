#include "GhostCellImmersedBoundaryNormalGradientBcStencil.h"

GhostCellImmersedBoundary::NormalGradientBcStencil::NormalGradientBcStencil(const Cell &cell,
                                                                            const ImmersedBoundaryObject &ibObj,
                                                                            const FiniteVolumeGrid2D &grid)
    :
      BcStencil(cell, ibObj, grid)
{
    initInterpolationCoeffs(ibObj, grid);
}

GhostCellImmersedBoundary::NormalGradientBcStencil::NormalGradientBcStencil(const Cell &cell,
                                                                            const ImmersedBoundaryObject &ibObj,
                                                                            const FiniteVolumeGrid2D &grid, const Vector2D &r)
    :
      BcStencil(cell, ibObj, grid, r)
{
    initInterpolationCoeffs(ibObj, grid);
}

void GhostCellImmersedBoundary::NormalGradientBcStencil::initInterpolationCoeffs(const ImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid)
{
    bool cellInStencil = false;

    for(const Cell& cell: cells_)
        if(cell_.get().id() == cell.id())
        {
            cellInStencil = true;
            break;
        }

    if(cellInStencil)
    {
        coeffs_ = StaticMatrix<1, 4>({bp_.y * nw_.x + bp_.x * nw_.y, nw_.x, nw_.y, 0.}) * Ainv_;
    }
    else
    {
        cells_.push_back(std::cref(cell_));

        Scalar l = length();
        coeffs_ = -StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Ainv_ / l;
        coeffs_.push_back(1. / l);
    }
}
