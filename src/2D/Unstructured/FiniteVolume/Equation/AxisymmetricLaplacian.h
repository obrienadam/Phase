#ifndef AXISYMMETRIC_LAPLACIAN_H
#define AXISYMMETRIC_LAPLACIAN_H

#include "Equation.h"

namespace axi
{
    Equation<Scalar> laplacian(Scalar gamma,
                               ScalarFiniteVolumeField &phi,
                               Scalar theta = 1.);

    Equation<Vector2D> vectorLaplacian(Scalar gamma,
                                       VectorFiniteVolumeField &phi,
                                       Scalar theta = 1.);
}

#endif
