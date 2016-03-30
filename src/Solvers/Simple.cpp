#include "Simple.h"
#include "Exception.h"

Simple::Simple(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(grid.addVectorField(input, "u")),
      h(grid.addVectorField(input, "h")),
      p(grid.addScalarField(input, "p")),
      pCorr(grid.addScalarField("pCorr")),
      rho(grid.addScalarField("rho")),
      mu(grid.addScalarField("mu")),
      m(grid.addScalarField("m")),
      d(grid.addScalarField("d")),
      uEqn_(grid, u),
      pCorrEqn_(grid, pCorr)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho"));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu"));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));
}

Scalar Simple::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    //solvePCorrEqn();
    return 0.;
}

//- Protected methods

Scalar Simple::solveUEqn(Scalar timeStep)
{
    interpolateFaces(u);
    uEqn_ = (rho*ddt(u, timeStep) + rho*div(u, u) == mu*laplacian(u) - grad(p));
    return uEqn_.solve();
}

Scalar Simple::solvePCorrEqn()
{
    pCorrEqn_ = (rho*d*laplacian(pCorr) == m);
    return pCorrEqn_.solve();
}

void Simple::computeD()
{
    const auto diag = uEqn_.matrix().diagonal();

    for(const Cell& cell: d.grid.cells)
        d[cell.id()] = cell.volume()/diag[cell.id()];
}

void Simple::rhieChowInterpolation()
{
    computeD();

    h = u;
}
