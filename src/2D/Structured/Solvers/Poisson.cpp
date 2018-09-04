#include "FiniteVolume/Discretization/Laplacian.h"

#include "Poisson.h"

Poisson::Poisson(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid)
    :
      Solver(input, grid),
      _phi(addField<Scalar>("phi")),
      _phiEqn("phiEqn", _phi)
{
    _gamma = input.caseInput().get<Scalar>("Properties.gamma", 1.);
}

Scalar Poisson::solve(Scalar timeStep)
{
    _phiEqn = (fv::lap(_gamma, _phi));
    return _phiEqn.solve();
}
