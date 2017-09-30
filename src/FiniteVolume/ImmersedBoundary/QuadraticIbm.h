#ifndef QUADRATIC_FACE_IBM_H
#define QUADRATIC_FACE_IBM_H

#include "FiniteVolumeEquation.h"
#include "ImmersedBoundary.h"

namespace qibm {

    class QuadraticStencil
    {
    public:

        QuadraticStencil(const Cell&, const ImmersedBoundaryObject& ibObj);

    private:

        std::vector<Ref<const Cell>> cells_;
        std::vector<Scalar> coeffs_;

    };

    Equation<Vector2D> div(const VectorFiniteVolumeField &phi, VectorFiniteVolumeField &u, const ImmersedBoundary &ib);

    Equation<Vector2D> laplacian(Scalar mu, VectorFiniteVolumeField &u, const ImmersedBoundary &ib);

    void computeFaceVelocities(VectorFiniteVolumeField &u, const ImmersedBoundary &ib);
}

#endif