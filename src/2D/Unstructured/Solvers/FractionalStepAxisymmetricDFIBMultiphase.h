#ifndef PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_MULTIPHASE_H
#define PHASE_FRACTIONAL_STEP_AXISYMMETRIC_DFIB_MULTIPHASE_H

#include "FractionalStepAxisymmetricDFIB.h"

class FractionalStepAxisymmetricDFIBMultiphase: public FractionalStepAxisymmetricDFIB
{
public:

    FractionalStepAxisymmetricDFIBMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

protected:

    virtual Scalar solveUEqn(Scalar timeStep) override;

    virtual Scalar solvePEqn(Scalar timeStep) override;

    virtual void correctVelocity(Scalar timeStep) override;

    virtual void computeIbForces(Scalar timeStep) override;

    Scalar rho1_, rho2_, mu1_, mu2_;

    ScalarFiniteVolumeField &gamma_, &rho_, &mu_;

    FiniteVolumeEquation<Scalar> gammaEqn_;
};

#endif
