#ifndef PHASE_AXISYMMETRIC_LAPLACIAN_H
#define PHASE_AXISYMMETRIC_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi {
FiniteVolumeEquation<Scalar> laplacian(Scalar gamma,
                                       ScalarFiniteVolumeField &phi);

FiniteVolumeEquation<Scalar> laplacian(const ScalarFiniteVolumeField &gamma,
                                       ScalarFiniteVolumeField &phi);

FiniteVolumeEquation<Vector2D>
laplacian(Scalar gamma, VectorFiniteVolumeField &u, Scalar theta = 1.);

FiniteVolumeEquation<Vector2D> laplacian(const ScalarFiniteVolumeField &mu,
                                         VectorFiniteVolumeField &u,
                                         Scalar theta = 1.);
} // namespace axi

#endif
