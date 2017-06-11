#ifndef FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H
#define FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H

#include "FractionalStepSimple.h"

class FractionalStepSimpleMultiphase: public FractionalStepSimple
{
public:
    FractionalStepSimpleMultiphase(const Input& input, const Communicator& comm, FiniteVolumeGrid2D& grid);

    void initialize();

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &rho, &mu, &gamma;
    VectorFiniteVolumeField &rhoU, &gradGamma;

private:

    Scalar solveGammaEqn(Scalar timeStep);

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_;

    Equation<Scalar> gammaEqn_;
};

#endif