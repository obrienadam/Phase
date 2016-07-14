#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "Equation.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> laplacian(ScalarFiniteVolumeField& field);
Equation<ScalarFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& coeff, ScalarFiniteVolumeField& field);
Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& coeff, VectorFiniteVolumeField& field);

}

#endif
