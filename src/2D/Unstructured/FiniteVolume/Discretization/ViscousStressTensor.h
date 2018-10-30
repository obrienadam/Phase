#ifndef PHASE_VISCOUS_STRESS_TENSOR_H
#define PHASE_VISCOUS_STRESS_TENSOR_H

#include "FiniteVolume/Equation/VectorFiniteVolumeEquation.h"

namespace fv {

VectorFvmEquation divTau(const ScalarFiniteVolumeField &mu, VectorFiniteVolumeField &u, Scalar theta = 1.);

}

#endif
