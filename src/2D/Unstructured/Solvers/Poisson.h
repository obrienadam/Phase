#ifndef PHASE_POISSON_H
#define PHASE_POISSON_H

#include "FiniteVolume/Equation/Equation.h"

#include "Solver.h"

class Poisson : public Solver
{
public:

    Poisson(const Input &input);

    Scalar solve(Scalar timeStep);

    Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
    { return std::numeric_limits<Scalar>::infinity(); }

    ScalarFiniteVolumeField &phi, &gamma;

protected:

    Equation<Scalar> phiEqn_;

};

#endif
