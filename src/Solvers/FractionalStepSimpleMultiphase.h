#ifndef FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H
#define FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H

#include "FractionalStepSimple.h"

class FractionalStepSimpleMultiphase: public FractionalStepSimple
{
public:
    FractionalStepSimpleMultiphase(const Input& input, const Communicator& comm, FiniteVolumeGrid2D& grid);

    ScalarFiniteVolumeField &rho, &mu, &gamma;
    VectorFiniteVolumeField &rhoU;

private:

    Equation<Scalar> gammaEqn_;
};

#endif