#ifndef FRACTIONAL_STEP_BOUSSINESQ_H
#define FRACTIONAL_STEP_BOUSSINESQ_H

#include "FractionalStep.h"

class FractionalStepBoussinesq : public FractionalStep
{
public:

    FractionalStepBoussinesq(const Input &input);

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &T;

protected:

    Equation<Scalar> TEqn_;

    Scalar solveUEqn(Scalar timeStep);

    Scalar solveTEqn(Scalar timeStep);

    Scalar alpha_, T0_, kappa_;

};


#endif
