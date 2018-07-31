#ifndef PHASE_FRACTIONAL_STEP_EULER_LAGRANGE_H
#define PHASE_FRACTIONAL_STEP_EULER_LAGRANGE_H

#include "FractionalStep.h"

class FractionalStepELIB : public FractionalStep
{
public:
    FractionalStepELIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    //virtual void correctVelocity(Scalar timeStep);
};


#endif //PHASE_FRACTIONALSTEPLAGRANGEEULER_H
