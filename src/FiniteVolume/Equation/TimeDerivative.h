#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include "Equation.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep);
Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep);
Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& rho, VectorFiniteVolumeField& field, Scalar timeStep);

}

#endif
