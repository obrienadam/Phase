#ifndef PHASE_POISSON_H
#define PHASE_POISSON_H

#include "Solver.h"

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

class Poisson: public Solver
{
public:

    Poisson(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep) override;

protected:

    ScalarField &_phi;

    FiniteVolumeEquation<Scalar> _phiEqn;

    Scalar _gamma;
};

#endif
