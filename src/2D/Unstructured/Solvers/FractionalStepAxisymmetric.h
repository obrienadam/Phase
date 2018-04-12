#ifndef PHASE_FRACTIONAL_STEP_AXISYMMETRIC_H
#define PHASE_FRACTIONAL_STEP_AXISYMMETRIC_H

#include "FractionalStep.h"

class FractionalStepAxisymmetric : public FractionalStep
{
public:

    FractionalStepAxisymmetric(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

protected:

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar maxDivergenceError();
};


#endif
