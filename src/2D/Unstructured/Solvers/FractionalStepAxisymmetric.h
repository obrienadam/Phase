#ifndef PHASE_FRACTIONAL_STEP_AXISYMMETRIC_H
#define PHASE_FRACTIONAL_STEP_AXISYMMETRIC_H

#include "FractionalStep.h"

class FractionalStepAxisymmetric : public FractionalStep
{
public:

    FractionalStepAxisymmetric(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    virtual Scalar maxCourantNumber(Scalar timeStep) const override;

protected:

    virtual Scalar solveUEqn(Scalar timeStep) override;

    virtual Scalar solvePEqn(Scalar timeStep) override;

    virtual void correctVelocity(Scalar timeStep) override;

    Scalar maxDivergenceError() override;
};


#endif
