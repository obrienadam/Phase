#include "Multiphase.h"
#include "Cicsam.h"
#include "Plic.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gamma", SparseMatrix::IncompleteLUT)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    setInitialConditions(input);

    //- Configuration
    interfaceAdvectionMethod_ = CICSAM;
    curvatureEvaluationMethod_ = CSF;

    if(interfaceAdvectionMethod_ == PLIC)
    {
        addGeometries("plicPolygons");
        addGeometries("fluxPolygons");
    }

    switch(curvatureEvaluationMethod_)
    {
    case CSF:
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, gamma, u, scalarFields_, vectorFields_));
        break;

    case HF:
        throw Exception("Multiphase", "Multiphase", "height functions not yet implemented.");
        break;

    }

    surfaceTensionForce_->compute();
    computeRho();
    computeMu();
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
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_->gammaTilde();

    for(const Cell& cell: rho.grid.activeCells())
    {
        size_t id = cell.id();
        rho[id] = (1. - gammaTilde[id])*rho1_ + gammaTilde[id]*rho2_;
    }

    //interpolateFaces(rho);
    harmonicInterpolateFaces(rho);
}

void Multiphase::computeMu()
{
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_->gammaTilde();

    for(const Cell& cell: mu.grid.activeCells())
    {
        size_t id = cell.id();
        mu[id] = (1. - gammaTilde[id])*mu1_ + gammaTilde[id]*mu2_;
    }

    harmonicInterpolateFaces(mu);
}

Scalar Multiphase::solveUEqn(Scalar timeStep)
{
    sg = fv::gravity(rho, g_);
    ft = surfaceTensionForce_->compute(); // surface tension forceb

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u)
             == fv::laplacian(mu, u) - fv::grad(p) + fv::source(ft) + fv::source(sg));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Multiphase::solveGammaEqn(Scalar timeStep)
{ 
    gamma.savePreviousTimeStep(timeStep, 1);
    interpolateFaces(gamma);

    switch(interfaceAdvectionMethod_)
    {
    case CICSAM:

        gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) == 0.);
        break;
    case PLIC:

        gammaEqn_ = (plic::div(u, gamma, timeStep, geometries()["plicPolygons"]) == 0.);
        break;
    }

    return gammaEqn_.solve();
}

void Multiphase::rhieChowInterpolation()
{
    Simple::rhieChowInterpolation();

    for(const Face& face: u.grid.interiorFaces())
    {
        const size_t id = face.id();
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar df = d.faces()[id];
        const Scalar rhof = rho.faces()[id];
        const Scalar rhoP = rho[lCell.id()];
        const Scalar rhoQ = rho[rCell.id()];

        Scalar kf;

        if(!u.grid.fluidCells().isInGroup(lCell))
            kf = surfaceTensionForce_->kappa()[rCell.id()];
        else if(!u.grid.fluidCells().isInGroup(rCell))
            kf = surfaceTensionForce_->kappa()[lCell.id()];
        else
            kf = surfaceTensionForce_->kappa().faces()[id];

        const Scalar gP = gamma[lCell.id()];
        const Scalar gQ = gamma[rCell.id()];

        u.faces()[id] += surfaceTensionForce_->sigma()*df*kf*(gQ - gP)*rc/dot(rc, rc) - rhof*(d[lCell.id()]*ft[lCell.id()]/rhoP + d[rCell.id()]*ft[rCell.id()]/rhoQ)/2.;
    }
}

//- External functions
