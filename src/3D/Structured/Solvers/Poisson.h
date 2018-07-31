#ifndef PHASE_POISSON_H
#define PHASE_POISSON_H

#include "Structured/FiniteVolume/Equation/ScalarFiniteVolumeEquation.h"

#include "Solver.h"

class Poisson: public Solver
{
public:

    Poisson(const Input &input, const std::shared_ptr<const StructuredGrid3D> &grid);

    Scalar solve(Scalar timeStep);

protected:

    ScalarField &_phi;

    ScalarFiniteVolumeEquation _phiEqn;

};

#endif
