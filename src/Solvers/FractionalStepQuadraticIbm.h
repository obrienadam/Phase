#ifndef FRACTIONAL_STEP_QUADRATIC_IBM_H
#define FRACTIONAL_STEP_QUADRATIC_IBM_H

#include "FractionalStep.h"

class FractionalStepQuadraticIbm: public FractionalStep
{
public:
    FractionalStepQuadraticIbm(const Input& input,
                   std::shared_ptr<FiniteVolumeGrid2D> &grid);

protected:
    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);
};


#endif
