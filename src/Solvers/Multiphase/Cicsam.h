#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace cicsam
{

Equation<Scalar> cn(const VectorFiniteVolumeField &u,
                    const VectorFiniteVolumeField& gradGamma,
                    const VectorFiniteVolumeField &m,
                    ScalarFiniteVolumeField &gamma,
                    Scalar timeStep);

}

#endif
