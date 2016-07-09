#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "Equation.h"

namespace cn
{

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field);

}

#endif
