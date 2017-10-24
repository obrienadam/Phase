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

    Type type() const
    { return GHOST_CELL; }

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField& u) const;

    Equation<Scalar> pressureBcs(Scalar rho, ScalarFiniteVolumeField& p) const;

    Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField& gamma, Scalar theta) const;

    const std::vector<GhostCellStencil>& stencils() const
    { return stencils_; }

protected:

    void constructStencils();

    std::vector<GhostCellStencil> stencils_;
};

#endif
