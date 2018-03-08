#ifndef FRACTIONAL_STEP_DIRECT_FORCING_H
#define FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"

class FractionalStepDirectForcing : public FractionalStep
{
public:
    using FractionalStep::FractionalStep;

protected:

    Scalar solveUEqn(Scalar timeStep);
};


#endif
