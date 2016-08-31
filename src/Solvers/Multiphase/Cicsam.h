#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace cicsam
{

enum Type{HC, UQ};

Scalar hc(Scalar gammaTilde, Scalar coD);
Scalar uq(Scalar gammaTilde, Scalar coD);

Equation<ScalarFiniteVolumeField> cn(const VectorFiniteVolumeField& u,
                                     const VectorFiniteVolumeField& gradGamma,
                                     ScalarFiniteVolumeField& gamma,
                                     Scalar timeStep,
                                     Type type = HC);

Equation<ScalarFiniteVolumeField> cn(const VectorFiniteVolumeField &u,
                                     const VectorFiniteVolumeField& gradGamma,
                                     const VectorFiniteVolumeField &m,
                                     ScalarFiniteVolumeField &gamma,
                                     Scalar timeStep);

}

#endif
