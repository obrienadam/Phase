#ifndef POISSON_H
#define POISSON_H

#include "Solver.h"
#include "Equation.h"

class Poisson : public Solver
{
public:

    Poisson(const FiniteVolumeGrid2D& grid, const Input& input);
    virtual Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField phi;

protected:

    Scalar gamma_;
    Equation<ScalarFiniteVolumeField> phiEqn_;
};

#endif
