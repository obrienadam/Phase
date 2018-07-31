#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H

#include "FractionalStepDFIB.h"

#include "Unstructured/FiniteVolume/Multiphase/CelesteImmersedBoundary.h"

class FractionalStepDirectForcingMultiphase: public FractionalStepDFIB
{
public:
    FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar solve(Scalar timeStep) override;

protected:

    Scalar solveGammaEqn(Scalar timeStep);

    Scalar solveUEqn(Scalar timeStep) override;

    Scalar solvePEqn(Scalar timeStep) override;

    void updateProperties(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_, capillaryTimeStep_;

    ScalarFiniteVolumeField &gamma_, &rho_, &mu_;

    VectorFiniteVolumeField &sg_, rhoU_;

    ScalarGradient &gradGamma_, &gradRho_;

    CelesteImmersedBoundary fst_;

    FiniteVolumeEquation<Scalar> gammaEqn_;
};

#endif
