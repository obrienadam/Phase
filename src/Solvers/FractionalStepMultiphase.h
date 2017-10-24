#ifndef FRACTIONAL_STEP_MULTIPHASE_H
#define FRACTIONAL_STEP_MULTIPHASE_H

#include "FractionalStep.h"
#include "Celeste.h"

class FractionalStepMultiphase : public FractionalStep
{
public:
    FractionalStepMultiphase(const Input &input,
                             std::shared_ptr<FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    virtual Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &rho, &mu, &gamma;
    ScalarGradient &gradGamma, &gradRho;
    VectorFiniteVolumeField &rhoU, &sg;
    Celeste &ft;

protected:

    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    void computeMomentumFlux(const ScalarFiniteVolumeField &beta, Scalar timeStep);

    void updateProperties(Scalar timeStep);

    //- Properties
    Scalar rho1_, rho2_, mu1_, mu2_, capillaryTimeStep_;

    //- Equations
    Equation<Scalar> gammaEqn_;
};

#endif