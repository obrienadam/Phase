#include "Multiphase.h"
#include "Cicsam.h"
#include "Plic.h"
#include "Celeste.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gamma", SparseMatrix::NoPreconditioner)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    setInitialConditions(input);

    //- Configuration
    interfaceAdvectionMethod_ = CICSAM;

    std::string tmp = input.caseInput().get<std::string>("Solver.surfaceTensionModel", "CSF");

    if(tmp == "CSF")
    {
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, gamma, u, scalarFields_, vectorFields_));
    }
    else if(tmp == "CELESTE")
    {
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new Celeste(input, gamma, u, scalarFields_, vectorFields_));
    }
    else
        throw Exception("Multiphase", "Multiphase", "unrecognized surface tension model \"" + tmp + "\".");

    //surfaceTensionForce_->compute();
    computeRho();
    computeMu();
}

Scalar Multiphase::solve(Scalar timeStep)
{
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
    ft = surfaceTensionForce_->compute(); // surface tension force. The order matters here since the smoothed gamma field is used for rho and mu
    computeRho();
    computeMu();
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u)
             == ab::laplacian(mu, u) - fv::grad(p) + fv::source(ft) + fv::source(sg));
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

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        u.faces()[id] += surfaceTensionForce_->sigma()*df*kf*(gQ - gP)*rc/dot(rc, rc) - rhof*(g*d[lCell.id()]*ft[lCell.id()]/rhoP + (1. - g)*d[rCell.id()]*ft[rCell.id()]/rhoQ)/2.;
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell& cellP = face.lCell();
        const Vector2D rf = face.centroid() - cellP.centroid();
        const Scalar df = d.faces()[face.id()];
        const Scalar kf = surfaceTensionForce_->kappa().faces()[face.id()];
        const Scalar gf = gamma.faces()[face.id()];
        const Scalar gP = gamma[cellP.id()];
        const Scalar rhoP = rho[cellP.id()];
        const Scalar rhof = rho.faces()[face.id()];

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u.faces()[face.id()] += surfaceTensionForce_->sigma()*df*kf*(gf - gP)*rf/dot(rf, rf)
                    - rhof*d[cellP.id()]*ft[cellP.id()]/rhoP;
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            break;

        case VectorFiniteVolumeField::OUTFLOW:
            u.faces()[face.id()] += dot(u.faces()[face.id()], face.outwardNorm(cellP.centroid())) > 0. ? surfaceTensionForce_->sigma()*df*kf*(gf - gP)*rf/dot(rf, rf)
                                                                                                         - rhof*d[cellP.id()]*ft[cellP.id()]/rhoP : Vector2D(0., 0.);
            break;

        default:
            throw Exception("Multiphase", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

//- External functions
