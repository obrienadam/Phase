#ifndef FRACTIONAL_STEP_MULTIPHASE_EULER_LAGRANGE_H
#define FRACTIONAL_STEP_MULTIPHASE_EULER_LAGRANGE_H

#include "FractionalStepEulerLagrange.h"

class FractionalStepMultiphaseEulerLagrange : public FractionalStepEulerLagrange
{
public:

    FractionalStepMultiphaseEulerLagrange(const Input &input,
                                          std::shared_ptr<FiniteVolumeGrid2D> &grid);
};


#endif
