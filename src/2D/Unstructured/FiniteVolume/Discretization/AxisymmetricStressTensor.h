#ifndef PHASE_AXISYMMETRIC_STRESS_TENSOR_H
#define PHASE_AXISYMMETRIC_STRESS_TENSOR_H

#include "FiniteVolume/Equation/VectorFiniteVolumeEquation.h"

namespace axi
{

VectorFvmEquation divSigma(Scalar mu, const ScalarFiniteVolumeField &p, VectorFiniteVolumeField &u, Scalar theta);

}

#endif
