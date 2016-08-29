#ifndef FRACTIONAL_STEP_MULTIPHASE_H
#define FRACTIONAL_STEP_MULTIPHASE_H

#include "FractionalStep.h"
#include "Celeste.h"

class FractionalStepMultiphase : public FractionalStep
{
public:

    FractionalStepMultiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);


    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &gradGamma, &ft;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void balancedForceInterpolation(Scalar timeStep);

    void computeRho();
    void computeMu();

    Scalar rho1_, rho2_, mu1_, mu2_;
    Celeste surfaceTensionForce_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;
};

#endif
