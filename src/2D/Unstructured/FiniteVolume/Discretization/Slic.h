#ifndef SLIC_H
#define SLIC_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace slic {

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep);

}

#endif
