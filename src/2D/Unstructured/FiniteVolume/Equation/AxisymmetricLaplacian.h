#ifndef AXISYMMETRIC_LAPLACIAN_H
#define AXISYMMETRIC_LAPLACIAN_H

#include "FiniteVolumeEquation.h"

namespace axi
{
    FiniteVolumeEquation<Scalar> laplacian(Scalar gamma,
                               ScalarFiniteVolumeField &phi,
                               Scalar theta = 1.);

    FiniteVolumeEquation<Vector2D> vectorLaplacian(Scalar gamma,
                                       VectorFiniteVolumeField &phi,
                                       Scalar theta = 1.);
}

#endif
