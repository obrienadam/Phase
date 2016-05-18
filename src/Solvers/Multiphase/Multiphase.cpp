#include "Multiphase.h"
#include "Cicsam.h"
#include "Plic.h"
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
    constructSmoothingKernels();

    //- Configuration
    interfaceAdvectionMethod_ = CICSAM;
    curvatureEvaluationMethod_ = CSF;
}

Scalar Multiphase::solve(Scalar timeStep, Scalar prevTimeStep)
{
    computeRho();
    computeMu();

    Piso::solve(timeStep, prevTimeStep);

    solveGammaEqn(timeStep, prevTimeStep);

    return 0.; // just to get rid of warning
}

//- Protected methods

void Multiphase::constructSmoothingKernels()
{
    CellSearch cs(grid_.activeCells());

    cellRangeSearch_.clear();
    cellRangeSearch_.resize(grid_.cells().size());

    for(const Cell &cell: grid_.activeCells())
        cellRangeSearch_[cell.id()] = cs.rangeSearch(Circle(cell.centroid(), kernelWidth_));
}

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

Scalar Multiphase::solveUEqn(Scalar timeStep, Scalar prevTimeStep)
{
    computeInterfaceNormals();
    computeCurvature();

    uEqn_ = (fv::ddt(rho, u, timeStep, prevTimeStep) + fv::div(rho*u, u)
             == fv::laplacian(mu, u) - fv::grad(p) + fv::source(sigma_*kappa*grad(gammaTilde)) + fv::source(rho*g_));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Multiphase::solveGammaEqn(Scalar timeStep, Scalar prevTimeStep)
{ 
    gamma.save(1);
    interpolateFaces(gamma);

    switch(interfaceAdvectionMethod_)
    {
    case CICSAM:

        gammaEqn_ = (fv::ddt(gamma, timeStep, prevTimeStep) + cicsam::div(u, gamma, timeStep) == 0.);
        break;
    case PLIC:

        gammaEqn_ = (fv::ddt(gamma, timeStep, prevTimeStep) + plic::div(u, gamma) == 0.);
        break;
    }

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
        vec = vec.unitVec();
        Scalar magSqr = vec.magSqr();

        if(!(fabs(magSqr - 1.) < 1e-2) || isnan(magSqr))
            vec.x = vec.y = 0.;
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
