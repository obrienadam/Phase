#ifndef PHASE_AXISYMMETRIC_CICSAM_H
#define PHASE_AXISYMMETRIC_CICSAM_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi
{
namespace cicsam
{

std::vector<Scalar> faceInterpolationWeights(const VectorFiniteVolumeField &u,
                                             const ScalarFiniteVolumeField &gamma,
                                             const VectorFiniteVolumeField &gradGamma,
                                             Scalar timeStep);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &gamma,
                                 const std::vector<Scalar> &faceInterpolationWeights,
                                 Scalar theta,
                                 const CellGroup &cells);

void computeMomentumFlux(Scalar rho1,
                         Scalar rho2,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &gamma,
                         const std::vector<Scalar> &faceInterpolationWeights,
                         VectorFiniteVolumeField &rhoU);
}
}

#endif
