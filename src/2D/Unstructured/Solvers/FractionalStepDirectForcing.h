#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"

class FractionalStepDirectForcing : public FractionalStep
{
public:
    FractionalStepDirectForcing(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    Scalar solve(Scalar timeStep);

    VectorFiniteVolumeField &fb;

protected:

    void solveExtEqns();

    void reconstructVelocity(Scalar timeStep);

    FiniteVolumeEquation<Vector2D> extEqn_;

    Scalar solveUEqn(Scalar timeStep);

};


#endif
