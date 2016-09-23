#ifndef FRACTIONAL_STEP_MULTIPHASE_H
#define FRACTIONAL_STEP_MULTIPHASE_H

#include "FractionalStep.h"
#include "Celeste.h"

class FractionalStepMultiphase : public FractionalStep
{
public:

    FractionalStepMultiphase(const Input& input, FiniteVolumeGrid2D& grid);

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;


    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &gradGamma, &ft;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);
    virtual void computeMassSource(Scalar timeStep);

    void computeRho();
    void computeMu();

    Scalar rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    Celeste surfaceTensionForce_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;
};

#endif
