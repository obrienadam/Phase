#ifndef PHASE_HYDRODYNAMIC_STRESS_H
#define PHASE_HYDRODYNAMIC_STRESS_H

#include "Types/Types.h"
#include "Geometry/Tensor2D.h"

class HydrodynamicStress
{
public:

    Scalar rho, p, mu;

    Tensor2D gradU;
};

#endif
