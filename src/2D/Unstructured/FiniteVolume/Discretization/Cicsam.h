#ifndef PHASE_CICSAM_H
#define PHASE_CICSAM_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace cicsam
{
Scalar hc(Scalar gammaDTilde, Scalar coD);

Scalar uq(Scalar gammaDTilde, Scalar coD);

std::vector<Scalar> faceInterpolationWeights(const VectorFiniteVolumeField &u,
                                             const ScalarFiniteVolumeField &gamma,
                                             const VectorFiniteVolumeField &gradGamma,
                                             Scalar timeStep);

void computeMomentumFlux(Scalar rho1,
                         Scalar rho2,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &gamma,
                         const std::vector<Scalar> &faceInterpolationWeights,
                         VectorFiniteVolumeField &rhoU);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &gamma,
                                 const std::vector<Scalar> &faceInterpolationWeights,
                                 Scalar theta,
                                 const CellGroup &cells);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &gamma,
                                 const std::vector<Scalar> &faceInterpolationWeights,
                                 Scalar theta);
}

#endif
