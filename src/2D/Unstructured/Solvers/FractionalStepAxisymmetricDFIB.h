#ifndef PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_H
#define PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_H

#include "FractionalStepAxisymmetric.h"

class DirectForcingImmersedBoundary;

class FractionalStepAxisymmetricDFIB: public FractionalStepAxisymmetric
{
public:

    FractionalStepAxisymmetricDFIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep) override;

protected:

    virtual Scalar solveUEqn(Scalar timeStep) override;

    virtual void computeIbForces(Scalar timeStep);

    std::shared_ptr<DirectForcingImmersedBoundary> ib_;

    VectorFiniteVolumeField &fib_;

    FiniteVolumeEquation<Vector2D> fibEqn_;
};

#endif
