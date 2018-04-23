#ifndef PHASE_CICSAM_H
#define PHASE_CICSAM_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace cicsam
{
    void beta(const VectorFiniteVolumeField &u,
              const VectorFiniteVolumeField &gradGamma,
              const ScalarFiniteVolumeField &gamma,
              Scalar timeStep,
              ScalarFiniteVolumeField &beta,
              Scalar k = 1.);

    void computeMomentumFlux(Scalar rho1,
                             Scalar rho2,
                             const VectorFiniteVolumeField &u,
                             const ScalarFiniteVolumeField &gamma,
                             const ScalarFiniteVolumeField &beta,
                             Scalar timeStep,
                             VectorFiniteVolumeField &rhoU);

    FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta);
}

#endif
