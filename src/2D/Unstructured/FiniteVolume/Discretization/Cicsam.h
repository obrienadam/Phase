#ifndef PHASE_CICSAM_H
#define PHASE_CICSAM_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace cicsam
{
std::vector<Scalar> beta(const VectorFiniteVolumeField &u,
                         const VectorFiniteVolumeField &gradGamma,
                         const ScalarFiniteVolumeField &gamma,
                         Scalar timeStep,
                         Scalar k = 1.);

void computeMomentumFlux(Scalar rho1,
                         Scalar rho2,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &gamma,
                         const std::vector<Scalar> &beta,
                         Scalar timeStep,
                         VectorFiniteVolumeField &rhoU);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 const std::vector<Scalar> &beta,
                                 ScalarFiniteVolumeField &gamma,
                                 Scalar theta,
                                 const CellGroup &cells);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 const std::vector<Scalar> &beta,
                                 ScalarFiniteVolumeField &gamma,
                                 Scalar theta);
}

#endif
