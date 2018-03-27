#ifndef HRIC_H
#define HRIC_H

#include "FiniteVolume/Equation/Equation.h"

namespace hric
{
    ScalarFiniteVolumeField beta(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 const ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep);

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField& beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta = 0.5);
}

#endif
