#ifndef FRACTIONAL_STEP_MULTIPHASE_H
#define FRACTIONAL_STEP_MULTIPHASE_H

#include "FractionalStep.h"
#include "Celeste.h"

class FractionalStepMultiphase : public FractionalStep
{
public:

    FractionalStepMultiphase(const Input& input, const Communicator &comm, FiniteVolumeGrid2D& grid);

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;


    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &gradGamma, &ft, &sg, &gradRho;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);
    virtual void correctVelocity(Scalar timeStep);

    void computeRho();
    void computeMu();

    Vector2D g_;

    Scalar rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    Celeste surfaceTensionForce_;

    Equation<Scalar> gammaEqn_;
};

#endif
