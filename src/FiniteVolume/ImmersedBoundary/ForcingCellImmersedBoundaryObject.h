#ifndef FORCING_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define FORCING_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "ForcingCellStencil.h"

class ForcingCellImmersedBoundaryObject: public ImmersedBoundaryObject
{
public:

    ForcingCellImmersedBoundaryObject(const std::string &name,
                                      Label id,
                                      FiniteVolumeGrid2D &grid);

    void update(Scalar timeStep);

    void updateCells();

    Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const;

protected:

    void constructStencils();

    CellZone pseudoFluidPoints_;

    std::vector<ForcingCellStencil> stencils_;
};

#endif
