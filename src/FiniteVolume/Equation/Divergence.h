#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "Equation.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField& u, ScalarFiniteVolumeField& field);
Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field);

}

#endif
