#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include "Equation.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, ScalarFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep);
Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep);
Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep);
Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& a, VectorFiniteVolumeField& field, Scalar timeStep, Scalar prevTimeStep);

}

#endif
