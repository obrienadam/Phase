#ifndef SLIC_H
#define SLIC_H

#include "Equation.h"

namespace slic {

    Equation<Scalar> div(const VectorFiniteVolumeField& u,
                         const VectorFiniteVolumeField& gradGamma,
                         ScalarFiniteVolumeField& gamma,
                         Scalar timeStep);

}

#endif

