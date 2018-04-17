#ifndef PHASE_STRUCTURED_FRACTIONAL_STEP_H
#define PHASE_STRUCTURED_FRACTIONAL_STEP_H

#include "StructuredSolver.h"

#include "FiniteDifference/FiniteDifferenceField.h"

class StructuredFractionalStep: public StructuredSolver
{
public:

    StructuredFractionalStep(const Input& input, const std::shared_ptr<StructuredGrid2D> &grid);

    Scalar solve(Scalar timeStep);

    FiniteDifferenceField<Vector2D> u;

    FiniteDifferenceField<Scalar> p;

protected:

    void solveUEqn(Scalar timeStep);

    void solvePEqn(Scalar timeStep);

    void correct(Scalar timeStep);

    Scalar _mu, _rho;
};


#endif
