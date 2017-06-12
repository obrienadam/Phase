#ifndef FRACTIONAL_STEP_MULTIPHASE_H
#define FRACTIONAL_STEP_MULTIPHASE_H

#include "FractionalStep.h"
#include "Celeste.h"

class FractionalStepMultiphase : public FractionalStep
{
public:

    FractionalStepMultiphase(const Input &input,
                             const Communicator &comm,
                             std::shared_ptr<FiniteVolumeGrid2D> &grid);

    virtual void initialize();

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;


    ScalarFiniteVolumeField &gamma, &rho, &mu;
    VectorFiniteVolumeField &gradGamma, &ft, &sg, &gradRho, &rhoU;

protected:

    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    void checkMassFluxConsistency(Scalar timeStep);

    void updateProperties(Scalar timeStep);

    Vector2D g_;

    Scalar cicsamBlending_, rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    Celeste surfaceTensionForce_;

    Equation<Scalar> gammaEqn_;
};

#endif
