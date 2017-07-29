#ifndef QUADRATIC_FACE_IBM_H
#define QUADRATIC_FACE_IBM_H

#include "FiniteVolumeEquation.h"
#include "ImmersedBoundary.h"

namespace qibm {
    Equation<Vector2D> div(const VectorFiniteVolumeField &phi, VectorFiniteVolumeField &u, const ImmersedBoundary &ib);

    Equation<Vector2D> laplacian(Scalar mu, VectorFiniteVolumeField &u, const ImmersedBoundary &ib);

    void computeFaceVelocities(VectorFiniteVolumeField &u, const ImmersedBoundary &ib);
}

#endif