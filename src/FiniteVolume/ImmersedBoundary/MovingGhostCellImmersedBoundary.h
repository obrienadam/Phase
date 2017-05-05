#ifndef MOVING_GHOST_CELL_IMMERSED_BOUNDARY_H
#define MOVING_GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundary.h"

namespace ib
{
    Equation<Vector2D>
    momentumEqn(const ImmersedBoundary& ib, const ScalarFiniteVolumeField &rho,
                const ScalarFiniteVolumeField &mu, const ScalarFiniteVolumeField &p, VectorFiniteVolumeField &u,
                Scalar timeStep);

    Equation<Scalar>
    pressureEqn(const ImmersedBoundary& ib, const ScalarFiniteVolumeField &rho,
                const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &p, Scalar timeStep);

    void correctVelocity(const ImmersedBoundary& ib, const ScalarFiniteVolumeField &rho,
                         ScalarFiniteVolumeField &p, VectorFiniteVolumeField& gradP, VectorFiniteVolumeField& u, Scalar timeStep);
}

#endif
