#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

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

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         const CellGroup &cells,
                         Scalar theta);

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta);
}

#endif
