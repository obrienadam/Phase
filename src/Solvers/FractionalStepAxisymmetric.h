#ifndef FRACTIONAL_STEP_AXISYMMETRIC_H
#define FRACTIONAL_STEP_AXISYMMETRIC_H

#include "FractionalStep.h"

class FractionalStepAxisymmetric : public FractionalStep
{
public:

    FractionalStepAxisymmetric(const Input &input,
                               std::shared_ptr<FiniteVolumeGrid2D> &grid);

protected:

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar maxDivergenceError();
};


#endif
