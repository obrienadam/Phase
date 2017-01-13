#include "Poisson.h"
#include "Laplacian.h"

Poisson::Poisson(const Input &input, FiniteVolumeGrid2D &grid)
    :
      Solver(input, grid),
      phi(addScalarField(input, "phi")),
      gamma(addScalarField("gamma")),
      phiEqn_(input, phi, "phiEqn")
{
    gamma.fill(input.caseInput().get<Scalar>("Properties.gamma", 1.));
    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (fv::laplacian(gamma, phi) + ib_.eqns(phi) == 0.);
    return phiEqn_.solve();
}
