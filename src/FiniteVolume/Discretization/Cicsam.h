#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

namespace cicsam
{
    Scalar hc(Scalar gammaDTilde, Scalar coD);

    Scalar uq(Scalar gammaDTilde, Scalar coD);

    ScalarFiniteVolumeField beta(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 const ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep,
                                 Scalar k = 1.);

    ScalarFiniteVolumeField massFluxes(Scalar rho1,
                                       Scalar rho2,
                                       const VectorFiniteVolumeField& u,
                                       const ScalarFiniteVolumeField& beta,
                                       const ScalarFiniteVolumeField& gamma,
                                       Scalar theta);

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField& beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta = 0.5);
}

#endif
