#include "Poisson.h"
#include "DiffusionTerm.h"

Poisson::Poisson(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      phi(grid.addScalarField(input, "phi")),
      phiEqn_(grid, phi)
{
    gamma_ = input.caseInput().get<Scalar>("Solver.gamma", 1.);
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (gamma_*laplacian(phi) == 0.);

    phiEqn_.solve();

    return phiEqn_.error();
}
