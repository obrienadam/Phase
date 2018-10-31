#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_MULTIPHASE_H

#include "FractionalStepDFIB.h"

class CelesteImmersedBoundary;

class FractionalStepDirectForcingMultiphase: public FractionalStepDFIB
{
public:
    FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar solve(Scalar timeStep) override;

protected:

    //- This class is used to communicate contact line info
    struct ContactLine
    {
        Point2D pt;

        Scalar beta, rho, rgh;

        Scalar gamma;

        Vector2D ncl, tcl;
    };

    Scalar solveGammaEqn(Scalar timeStep);

    Scalar solveUEqn(Scalar timeStep) override;

    Scalar solvePEqn(Scalar timeStep) override;

    void updateProperties(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    void computeIbForces(Scalar timeStep);

    void computeFieldExtenstions(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_, capillaryTimeStep_;

    ScalarFiniteVolumeField &gamma_, &rho_, &mu_, &gammaSrc_, &divU_;

    VectorFiniteVolumeField &sg_, &us_;

    ScalarGradient &gradGamma_, &gradRho_;

    std::shared_ptr<CelesteImmersedBoundary> fst_;

    FiniteVolumeEquation<Scalar> gammaEqn_;

    std::vector<ContactLine> contactLines_;
};

#endif
