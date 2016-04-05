#include "Poisson.h"
#include <stdio.h>

Poisson::Poisson(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      phi(grid.addScalarField(input, "phi")),
      gamma(grid.addScalarField("gamma")),
      phiEqn_(phi)
{
    gamma.fill(input.caseInput().get<Scalar>("Properties.gamma", 1.));
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (laplacian(gamma, phi) == 0.);
    phiEqn_.solve();
    return phiEqn_.error();
}
