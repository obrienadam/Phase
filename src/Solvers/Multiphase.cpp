#include "Multiphase.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Simple(grid, input),
      gamma(grid.addScalarField(input, "gamma"))
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");

    input.setInitialConditions(grid);

    computeRho();
    computeMu();
}

//- Private methods

void Multiphase::computeRho()
{
    for(const Cell& cell: rho.grid.cells)
    {
        size_t id = cell.id();
        rho[id] = rho1_*(1. - gamma[id]) + rho2_*gamma[id];
    }

    interpolateFaces(rho);
}

void Multiphase::computeMu()
{
    for(const Cell& cell: mu.grid.cells)
    {
        size_t id = cell.id();
        mu[id] = mu1_*(1. - gamma[id]) + mu2_*gamma[id];
    }

    interpolateFaces(mu);
}
