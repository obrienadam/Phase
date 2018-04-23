#ifndef PHASE_FRACTIONAL_STEP_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_MULTIPHASE_H

#include "FiniteVolume/Multiphase/Celeste.h"

#include "FractionalStep.h"

class FractionalStepMultiphase : public FractionalStep
{
public:
    FractionalStepMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    virtual Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &rho, &mu, &gamma, &beta;

    ScalarGradient &gradGamma, &gradRho;

    VectorFiniteVolumeField &rhoU, &sg;

    Celeste &ft;

protected:

    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    virtual void updateProperties(Scalar timeStep);

    //- Properties
    Scalar rho1_, rho2_, mu1_, mu2_, capillaryTimeStep_;

    //- Equations
    FiniteVolumeEquation<Scalar> gammaEqn_;
};

#endif