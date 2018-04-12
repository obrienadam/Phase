#ifndef PHASE_POISSON_H
#define PHASE_POISSON_H

#include "FiniteVolume/Equation/Equation.h"

#include "Solver.h"

class Poisson : public Solver
{
public:

    Poisson(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar solve(Scalar timeStep);

    Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
    { return std::numeric_limits<Scalar>::infinity(); }

    ScalarFiniteVolumeField &phi;

protected:

    Equation<Scalar> phiEqn_;

    CellGroup &solid_;

    Scalar gamma_;
};

#endif
