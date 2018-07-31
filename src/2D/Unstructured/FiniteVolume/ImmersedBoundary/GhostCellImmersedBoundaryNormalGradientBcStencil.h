#ifndef PHASE_GHOST_CELL_IMMERSED_BOUNDARY_NORMAL_GRADIENT_BC_STENCIL_H
#define PHASE_GHOST_CELL_IMMERSED_BOUNDARY_NORMAL_GRADIENT_BC_STENCIL_H

#include "GhostCellImmersedBoundaryBcStencil.h"

class GhostCellImmersedBoundary::NormalGradientBcStencil: public GhostCellImmersedBoundary::BcStencil
{
public:

    NormalGradientBcStencil(const Cell &cell,
                            const ImmersedBoundaryObject &ibObj,
                            const FiniteVolumeGrid2D &grid);

    NormalGradientBcStencil(const Cell &cell,
                            const ImmersedBoundaryObject &ibObj,
                            const FiniteVolumeGrid2D &grid,
                            const Vector2D &r);

protected:

    void initInterpolationCoeffs(const ImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid);
};

#endif
