#ifndef PHASE_AXISYMMETRIC_LAPLACIAN_H
#define PHASE_AXISYMMETRIC_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi
{
FiniteVolumeEquation<Scalar> laplacian(Scalar gamma,
                                       ScalarFiniteVolumeField &phi);

FiniteVolumeEquation<Vector2D> laplacian(Scalar gamma,
                                         VectorFiniteVolumeField &u,
                                         Scalar theta = 1.);
}

#endif
