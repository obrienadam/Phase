#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H

#include "Equation.h"

namespace ab
{

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& coeff, VectorFiniteVolumeField& field);

}

#endif
