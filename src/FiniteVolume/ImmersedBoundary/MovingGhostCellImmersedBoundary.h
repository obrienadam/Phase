#ifndef MOVING_GHOST_CELL_IMMERSED_BOUNDARY_H
#define MOVING_GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace ib
{
    Equation<Vector2D>
    ddt(const std::vector<Ref<const ImmersedBoundaryObject> > &ibObjs, const ScalarFiniteVolumeField &rho,
        VectorFiniteVolumeField &u, Scalar timeStep);

    Equation<Vector2D>
    div(const std::vector<Ref<const ImmersedBoundaryObject> > &ibObjs, const ScalarFiniteVolumeField &rho,
        VectorFiniteVolumeField &u, Scalar timeStep);
}

#endif