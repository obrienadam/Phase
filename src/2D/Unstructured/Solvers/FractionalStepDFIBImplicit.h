#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_IMPLICIT_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_IMPLICIT_H

#include "FractionalStepDFIB.h"

class FractionalStepDFIBImplicit : public FractionalStepDFIB
{
public:
    FractionalStepDFIBImplicit(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

protected:

    virtual Scalar solveUEqn(Scalar timeStep) override;

};


#endif
