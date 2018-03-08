#ifndef FRACTIONAL_STEP_EULER_LAGRANGE_H
#define FRACTIONAL_STEP_EULER_LAGRANGE_H

#include "FractionalStep.h"

class FractionalStepEulerLagrange : public FractionalStep
{
public:
    FractionalStepEulerLagrange(const Input &input);

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    //virtual void correctVelocity(Scalar timeStep);
};


#endif //PHASE_FRACTIONALSTEPLAGRANGEEULER_H
