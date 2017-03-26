#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

namespace cicsam
{

void interpolateFaces(const VectorFiniteVolumeField &u,
                      const VectorFiniteVolumeField &gradGamma,
                      ScalarFiniteVolumeField &gamma,
                      Scalar timeStep,
                      Scalar k = 1.);

Equation<Scalar> div(const VectorFiniteVolumeField &u,
                     ScalarFiniteVolumeField &gamma);

Equation<Scalar> div(const VectorFiniteVolumeField &u,
                     const VectorFiniteVolumeField &gradGamma,
                     ScalarFiniteVolumeField &gamma,
                     Scalar timeStep,
                     Scalar k = 1.);
}

#endif
