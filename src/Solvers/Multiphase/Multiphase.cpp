#include "Multiphase.h"
#include "Cicsam.h"
#include "Plic.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gamma")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    setInitialConditions(input);

    computeRho();
    computeMu();

    //- Configuration
    interfaceAdvectionMethod_ = PLIC;
    curvatureEvaluationMethod_ = CSF;

    if(interfaceAdvectionMethod_ == PLIC)
    {
        addGeometries("plicPolygons");
        addGeometries("fluxPolygons");
    }

    switch(curvatureEvaluationMethod_)
    {
    case CSF:

        surfaceTensionForce_ = std::unique_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, gamma));
        break;
    case HF:

        surfaceTensionForce_ = std::unique_ptr<SurfaceTensionForce>(new HeightFunction(input, gamma));
        interfaceAdvectionMethod_ = PLIC; // Only plic is guaranteed to work with height function methods
        addGeometries("hfPolygons");
        break;
    }
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

void Multiphase::computeRho()
{
//    for(const Cell& cell: rho.grid.activeCells())
//    {
//        size_t id = cell.id();
//        rho[id] = rho1_*(1. - gamma[id]) + rho2_*gamma[id];
//    }

//    harmonicInterpolateFaces(rho);

    for(const Cell& cell: rho.grid.activeCells())
    {
        size_t id = cell.id();
        rho[id] = (1. - gamma[id])*rho1_ + gamma[id]*rho2_;
    }

    interpolateFaces(rho);
}

void Multiphase::computeMu()
{
//    for(const Cell& cell: mu.grid.activeCells())
//    {
//        size_t id = cell.id();
//        mu[id] = rho[id]/((1. - gamma[id])*rho1_/mu1_ + gamma[id]*rho2_/mu2_);
//    }

//    interpolateFaces(mu);

    for(const Cell& cell: mu.grid.activeCells())
    {
        size_t id = cell.id();
        mu[id] = (1. - gamma[id])*mu1_ + gamma[id]*mu2_;
    }

    interpolateFaces(mu);
}

Scalar Multiphase::solveUEqn(Scalar timeStep, Scalar prevTimeStep)
{
    ft = surfaceTensionForce_->compute(); // surface tension force

    uEqn_ = (fv::ddt(rho, u, timeStep, prevTimeStep) + fv::div(rho*u, u)
             == fv::laplacian(mu, u) - fv::grad(p) + fv::source(ft) + fv::source(rho*g_));
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

        gammaEqn_ = (plic::div(u, gamma, timeStep, geometries()["plicPolygons"]) == 0.);
        break;
    }

    Scalar error = gammaEqn_.solve();

    if(isnan(error))
        throw Exception("Multiphase", "solveGammaEqn", "a nan value was detected.");

    return error;
}

//- External functions
