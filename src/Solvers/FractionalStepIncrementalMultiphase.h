#ifndef FRACTIONAL_STEP_INCREMENTAL_MULTIPHASE_H
#define FRACTIONAL_STEP_INCREMENTAL_MULTIPHASE_H

#include "FractionalStepIncremental.h"
#include "Celeste.h"

class FractionalStepIncrementalMultiphase : public FractionalStepIncremental
{
public:

    FractionalStepIncrementalMultiphase(const Input &input,
                             std::shared_ptr<FiniteVolumeGrid2D> &grid);

    virtual void initialize();

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;


    ScalarFiniteVolumeField &gamma, &rho, &mu;
    ScalarGradient &gradGamma, &gradRho;
    VectorFiniteVolumeField &sg, &rhoU;
    Celeste &ft;

protected:

    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    void updateProperties(Scalar timeStep);

    Vector2D g_;
    Scalar cicsamBlending_, rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    Equation<Scalar> gammaEqn_;
};

#endif
