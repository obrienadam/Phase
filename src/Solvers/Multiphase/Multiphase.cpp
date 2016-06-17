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

    computeRho();
    computeMu();

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
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, gamma, u, scalarFields_));
        break;

    case HF:
        throw Exception("Multiphase", "Multiphase", "height functions not yet implemented.");
        break;

    }
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
        rho[id] = (1. - gamma[id])*rho1_ + gamma[id]*rho2_;
    }

    //interpolateFaces(rho);
    harmonicInterpolateFaces(rho);
}

void Multiphase::computeMu()
{
    for(const Cell& cell: mu.grid.activeCells())
    {
        size_t id = cell.id();
        mu[id] = (1. - gamma[id])*mu1_ + gamma[id]*mu2_;
    }

    interpolateFaces(mu);
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

        gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, surfaceTensionForce_->n(), gamma, timeStep) == 0.);
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

        const Scalar kf = surfaceTensionForce_->kappa().faces()[id];
        const Scalar kP = surfaceTensionForce_->kappa()[lCell.id()];
        const Scalar kQ = surfaceTensionForce_->kappa()[rCell.id()];

        const Scalar gP = surfaceTensionForce_->gammaTilde()[lCell.id()];
        const Scalar gQ = surfaceTensionForce_->gammaTilde()[rCell.id()];

        const Vector2D& dgP = surfaceTensionForce_->gradGamma()[lCell.id()];
        const Vector2D& dgQ = surfaceTensionForce_->gradGamma()[rCell.id()];

        u.faces()[id] += df*surfaceTensionForce_->sigma()*(kf*(gQ - gP)*rc/dot(rc, rc) - rhof*(kP*dgP/rhoP + kQ*dgQ/rhoQ)/2.);
    }
}

//- External functions
