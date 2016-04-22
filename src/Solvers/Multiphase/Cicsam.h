#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

namespace hc
{

Scalar betaFace(Scalar gammaD, Scalar gammaA, Scalar gammaU, Scalar coD);

}

namespace sc
{

}

namespace uq
{

}

namespace cicsam
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep);

}

#endif
