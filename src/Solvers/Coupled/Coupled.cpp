#include "Coupled.h"

Coupled::Coupled(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      p(addScalarField(input, "p")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      eqn_(rho, mu, u, p)

{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho"));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu"));
}

Scalar Coupled::solve(Scalar timeStep)
{
    return eqn_.solve(timeStep);
}

Scalar Coupled::computeMaxTimeStep(Scalar maxCo) const
{
    Scalar maxTimeStep = std::numeric_limits<Scalar>::infinity();

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaX = nb.rCellVec().mag();
            Scalar magU = u.faces()[nb.face().id()].mag();

            maxTimeStep = std::min(maxTimeStep, maxCo*deltaX/magU);
        }

    return maxTimeStep;
}
