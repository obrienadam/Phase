#ifndef FRACTIONAL_STEP_AXISYMMETRIC_H
#define FRACTIONAL_STEP_AXISYMMETRIC_H

#include "FractionalStep.h"

class FractionalStepAxisymmetric : public FractionalStep
{
public:

    FractionalStepAxisymmetric(const Input &input);

protected:

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar maxDivergenceError();
};


#endif
