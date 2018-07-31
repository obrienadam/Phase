#ifndef PHASE_GHOST_CELL_IMMERSED_BOUNDARY_H
#define PHASE_GHOST_CELL_IMMERSED_BOUNDARY_H

#include "ImmersedBoundary.h"

class GhostCellImmersedBoundary: public ImmersedBoundary
{
public:

    class BcStencil;

    class FixedBcStencil;

    class NormalGradientBcStencil;

    GhostCellImmersedBoundary(const Input &input,
                              const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                              const std::shared_ptr<CellGroup> &domainCells);

    void updateCells() override;

    FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &phi) const;

    FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &u) const;

    FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    void computeForce(Scalar rho,
                      Scalar mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

protected:

};

#endif
