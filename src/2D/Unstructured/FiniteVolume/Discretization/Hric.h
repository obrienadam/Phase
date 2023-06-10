#ifndef HRIC_H
#define HRIC_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace hric {
ScalarFiniteVolumeField beta(const VectorFiniteVolumeField &u,
                             const VectorFiniteVolumeField &gradGamma,
                             const ScalarFiniteVolumeField &gamma,
                             Scalar timeStep);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 const ScalarFiniteVolumeField &beta,
                                 ScalarFiniteVolumeField &gamma,
                                 Scalar theta = 0.5);
} // namespace hric

#endif
