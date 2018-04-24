#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H

#include "FractionalStepMultiphase.h"

class FractionalStepDirectForcingMultiphase: public FractionalStepMultiphase
{
public:
    FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    Scalar solve(Scalar timeStep);

    VectorFiniteVolumeField &fb;

protected:

    Scalar solveUEqn(Scalar timeStep);
};

#endif
