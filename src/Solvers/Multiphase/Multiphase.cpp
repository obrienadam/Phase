#include "Multiphase.h"
#include "Cicsam.h"
#include "CellSearch.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      gammaTilde(addScalarField("gammaTilde")),
      kappa(addScalarField("kappa")),
      n(addVectorField("n")),
      gammaEqn_(gamma, "gamma")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");

    setInitialConditions(input);

    computeRho();
    computeMu();

    //- Construct the range search
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");
    cellRangeSearch_ = rangeSearch(grid, kernelWidth_);
}

Scalar Multiphase::solve(Scalar timeStep)
{
    computeRho();
    computeMu();

    Piso::solve(timeStep);

    solveGammaEqn(timeStep);

    return 0.; // just to get rid of warning
}

//- Protected methods

void Multiphase::computeRho()
{
    for(const Cell& cell: rho.grid.activeCells())
    {
        size_t id = cell.id();
        rho[id] = rho1_*(1. - gamma[id]) + rho2_*gamma[id];
    }

    harmonicInterpolateFaces(rho);
}

void Multiphase::computeMu()
{
    for(const Cell& cell: mu.grid.activeCells())
    {
        size_t id = cell.id();
        mu[id] = rho[id]/((1. - gamma[id])*rho1_/mu1_ + gamma[id]*rho2_/mu2_);
    }

    interpolateFaces(mu);
}

Scalar Multiphase::solveUEqn(Scalar timeStep)
{
    computeInterfaceNormals();
    computeCurvature();

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) == fv::laplacian(mu, u) - fv::grad(p) + fv::source(sigma_*kappa*grad(gammaTilde)) + fv::source(rho*g_));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Multiphase::solveGammaEqn(Scalar timeStep)
{
    interpolateFaces(gamma);
    gamma.save();
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, gamma, timeStep) == 0.);

    Scalar error = gammaEqn_.solve();

    if(isnan(error))
        throw Exception("Multiphase", "solveGammaEqn", "a nan value was detected.");

    return error;
}

void Multiphase::computeInterfaceNormals()
{
    gammaTilde = smooth(gamma, cellRangeSearch_, kernelWidth_);
    interpolateFaces(gammaTilde);
    n = grad(gammaTilde);

    for(Vector2D &vec: n)
    {
        if(vec.mag() > 1e-5)
            vec = vec.unitVec();
        else
            vec = Vector2D(0., 0.);
    }

    interpolateFaces(n);
}

void Multiphase::computeCurvature()
{
    for(const Cell &cell: kappa.grid.fluidCells())
    {
        Scalar &k = kappa[cell.id()] = 0.;

        for(const InteriorLink &nb: cell.neighbours())
            k -= dot(n.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            k -= dot(n.faces()[bd.face().id()], bd.outwardNorm());

        k /= cell.volume();
    }
}

//- External functions
