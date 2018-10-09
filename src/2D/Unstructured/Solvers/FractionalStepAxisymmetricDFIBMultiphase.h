#ifndef PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_MULTIPHASE_H

#include "FiniteVolume/Multiphase/CelesteAxisymmetricImmersedBoundary.h"

#include "FractionalStepAxisymmetricDFIB.h"

class FractionalStepAxisymmetricDFIBMultiphase: public FractionalStepAxisymmetricDFIB
{
public:

    FractionalStepAxisymmetricDFIBMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    virtual void initialize() override;

    virtual Scalar solve(Scalar timeStep) override;

protected:

    virtual Scalar solveGammaEqn(Scalar timeStep);

    void updateProperties(Scalar timeStep);

    virtual Scalar solveUEqn(Scalar timeStep) override;

    virtual Scalar solvePEqn(Scalar timeStep) override;

    virtual void correctVelocity(Scalar timeStep) override;

    virtual void computeIbForces(Scalar timeStep) override;

    virtual void computeFieldExtenstions(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_;

    ScalarFiniteVolumeField &gamma_, &gammaSrc_, &rho_, &mu_;

    VectorFiniteVolumeField &rhoU_, &sg_;

    ScalarGradient &gradGamma_, &gradRho_;

    CelesteAxisymmetricImmersedBoundary fst_;

    FiniteVolumeEquation<Scalar> gammaEqn_;
};

#endif
