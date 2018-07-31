#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"

class FractionalStepDFIB : public FractionalStep
{
public:
    FractionalStepDFIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep) override;

protected:

    virtual void solveExtEqns();

    virtual Scalar solveUEqn(Scalar timeStep) override;

    VectorFiniteVolumeField &fb_;

    FiniteVolumeEquation<Vector2D> extEqn_;

    std::shared_ptr<DirectForcingImmersedBoundary> ib_;
};


#endif
