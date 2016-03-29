#include "Simple.h"
#include "Exception.h"

Simple::Simple(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(grid.addVectorField(input, "u")),
      p(grid.addScalarField(input, "p")),
      pCorr(grid.addScalarField("pCorr")),
      rho(grid.addScalarField("rho")),
      mu(grid.addScalarField("mu")),
      m(grid.addScalarField("m")),
      uEqn_(grid, u),
      pCorrEqn_(grid, pCorr)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho"));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu"));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));
}

Scalar Simple::solve(Scalar timeStep)
{
    solveUEqn();
    //solvePCorrEqn();
    return 0.;
}

//- Protected methods

Scalar Simple::solveUEqn()
{
    uEqn_ = (mu*laplacian(u) == 0.);
    return uEqn_.solve();
}

Scalar Simple::solvePCorrEqn()
{
    pCorrEqn_ = (rho*laplacian(pCorr) == m);
    return pCorrEqn_.solve();
}
