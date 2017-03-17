#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

namespace cicsam
{

void interpolateFaces(const VectorFiniteVolumeField &u,
                      const VectorFiniteVolumeField &gradGamma,
                      const VectorFiniteVolumeField &m,
                      ScalarFiniteVolumeField &gamma,
                      Scalar timeStep, Scalar k = 1.);

Equation<Scalar> cn(const VectorFiniteVolumeField &u,
                    const VectorFiniteVolumeField &gradGamma,
                    const VectorFiniteVolumeField &m,
                    ScalarFiniteVolumeField &gamma,
                    Scalar timeStep,
                    Scalar alpha = 0.5, Scalar k = 1.);

}

#endif
