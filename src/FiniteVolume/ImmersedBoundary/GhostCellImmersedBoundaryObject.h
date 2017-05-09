#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "GhostCellStencil.h"

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

protected:

    void constructStencils();

    std::vector<GhostCellStencil> stencils_;
};

#endif
