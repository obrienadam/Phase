#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include "Equation.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, ScalarFiniteVolumeField& field, Scalar timeStep);
Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep);
Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep);
Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, VectorFiniteVolumeField& field, Scalar timeStep);

}

#endif
