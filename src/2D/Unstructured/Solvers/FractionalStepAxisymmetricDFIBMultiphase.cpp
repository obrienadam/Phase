#include "FractionalStepAxisymmetricDFIBMultiphase.h"

FractionalStepAxisymmetricDFIBMultiphase::FractionalStepAxisymmetricDFIBMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepAxisymmetricDFIB(input, grid),
      gamma_(*addField<Scalar>(input, "gamma", fluid_)),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      gammaEqn_(input, gamma_, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", FractionalStep::rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", FractionalStep::rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", FractionalStep::mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", FractionalStep::mu_);
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solveUEqn(Scalar timeStep)
{

}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solvePEqn(Scalar timeStep)
{

}

void FractionalStepAxisymmetricDFIBMultiphase::correctVelocity(Scalar timeStep)
{

}

void FractionalStepAxisymmetricDFIBMultiphase::computeIbForces(Scalar timeStep)
{

}
