#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "GhostCellStencil.h"
#include "Celeste.h"

class GhostCellImmersedBoundaryObject: public ImmersedBoundaryObject
{
public:

    GhostCellImmersedBoundaryObject(const std::string &name,
                                    Label id,
                                    FiniteVolumeGrid2D &grid);

    void update(Scalar timeStep);

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const;

    Equation<Vector2D> solidVelocity(VectorFiniteVolumeField& u) const;

    Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField& gamma, Scalar theta) const;

    Vector2D normalForce() const
    { return normalForce_; }

    Vector2D shearForce() const
    { return shearForce_; }

    void computeNormalForce(const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& p);

    void computeShearForce(const ScalarFiniteVolumeField &mu, const VectorFiniteVolumeField& u);

    const std::vector<GhostCellStencil>& stencils() const
    { return stencils_; }

protected:

    void constructStencils();

    std::vector<GhostCellStencil> stencils_;

    Vector2D normalForce_, shearForce_;
};

#endif
