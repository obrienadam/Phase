#ifndef POISSON_H
#define POISSON_H

#include "Solver.h"
#include "Communicator.h"
#include "Equation.h"

class Poisson : public Solver
{
public:

    Poisson(const Input &input,
            std::shared_ptr<FiniteVolumeGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
    { return std::numeric_limits<Scalar>::infinity(); }

    ScalarFiniteVolumeField &phi, &gamma;

protected:

    Equation<Scalar> phiEqn_;

};

#endif
