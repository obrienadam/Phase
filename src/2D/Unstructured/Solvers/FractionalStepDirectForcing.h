#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"

class FractionalStepDirectForcing : public FractionalStep
{
public:
    FractionalStepDirectForcing(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &divU;

    VectorFiniteVolumeField &fb;

protected:

    Equation<Scalar> pExtEqn_;

    Equation<Vector2D> uExtEqn_;

    void solveExtEqns();

    Scalar solveUEqn(Scalar timeStep);

};


#endif
