#ifndef PLIC_H
#define PLIC_H

#include "Equation.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field);

}

#endif
