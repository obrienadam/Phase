#ifndef FRACTIONAL_STEP_DIRECT_FORCING_H
#define FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"

class FractionalStepDirectForcing : public FractionalStep
{
public:
    FractionalStepDirectForcing(const Input &input);

    VectorFiniteVolumeField &fb;

protected:

    Scalar solveUEqn(Scalar timeStep);
};


#endif
