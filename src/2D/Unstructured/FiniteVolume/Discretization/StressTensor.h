#ifndef PHASE_STRESS_TENSOR
#define PHASE_STRESS_TENSOR

#include "FiniteVolume/Equation/VectorFiniteVolumeEquation.h"

namespace fv {

VectorFvmEquation divSigma(Scalar mu, const ScalarFiniteVolumeField &p,
                           VectorFiniteVolumeField &u, Scalar theta);

VectorFvmEquation divSigma(const ScalarFiniteVolumeField &mu,
                           const ScalarFiniteVolumeField &p,
                           VectorFiniteVolumeField &u, Scalar theta);

} // namespace fv

#endif
