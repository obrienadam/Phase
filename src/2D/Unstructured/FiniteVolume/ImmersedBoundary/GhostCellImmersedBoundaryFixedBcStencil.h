#ifndef PHASE_GHOST_CELL_IMMERSED_BOUNDARY_FIXED_BC_STENCIL_H
#define PHASE_GHOST_CELL_IMMERSED_BOUNDARY_FIXED_BC_STENCIL_H

#include "GhostCellImmersedBoundaryBcStencil.h"

class GhostCellImmersedBoundary::FixedBcStencil: public GhostCellImmersedBoundary::BcStencil
{
public:

    FixedBcStencil(const Cell &cell,
                   const ImmersedBoundaryObject &ibObj,
                   const FiniteVolumeGrid2D &grid);

    FixedBcStencil(const Cell &cell,
                   const ImmersedBoundaryObject &ibObj,
                   const FiniteVolumeGrid2D &grid,
                   const Vector2D &r);

protected:

    void initInterpolationCoeffs(const ImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid);
};

#endif
