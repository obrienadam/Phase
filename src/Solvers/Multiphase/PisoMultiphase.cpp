#include "PisoMultiphase.h"
#include "Cicsam.h"
#include "Plic.h"
#include "Celeste.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
#include "FaceInterpolation.h"

PisoMultiphase::PisoMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
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

    std::string tmp = input.caseInput().get<std::string>("Solver.surfaceTensionModel");

    if(tmp == "CSF")
    {
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, *this));
    }
    else if(tmp == "CELESTE")
    {
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new Celeste(input, *this));
    }
    else
        throw Exception("PisoMultiphase", "PisoMultiphase", "unrecognized surface tension model \"" + tmp + "\".");

    //surfaceTensionForce_->compute();
    computeRho();
    computeMu();
}

Scalar PisoMultiphase::solve(Scalar timeStep)
{
    Piso::solve(timeStep);
    solveGammaEqn(timeStep);

    return 0.; // just to get rid of warning
}

//- Protected methods

void PisoMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_->gammaTilde();

    for(const Cell& cell: rho.grid.activeCells())
    {
        rho(cell) = (1. - gammaTilde(cell))*rho1_ + gammaTilde(cell)*rho2_;
    }

    interpolateFaces(fv::INVERSE_VOLUME, rho);
}

void PisoMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_->gammaTilde();

    for(const Cell& cell: mu.grid.activeCells())
    {
        mu(cell) = (1. - gammaTilde(cell))*mu1_ + gammaTilde(cell)*mu2_;
        //mu[id] = rho[id]/((1. - gammaTilde[id])*rho1_/mu1_ + gammaTilde[id]*rho2_/mu2_);
    }

    interpolateFaces(fv::INVERSE_VOLUME, mu);
}

Scalar PisoMultiphase::solveUEqn(Scalar timeStep)
{
    ft = surfaceTensionForce_->compute(); // surface tension force. The order matters here since the smoothed gamma field is used for rho and mu
    computeRho();
    computeMu();
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_.eqns(u)
             == fv::laplacian(mu, u) - fv::source(gradP) + fv::source(ft) + fv::source(sg));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar PisoMultiphase::solveGammaEqn(Scalar timeStep)
{ 
    gamma.savePreviousTimeStep(timeStep, 1);
    interpolateFaces(fv::INVERSE_VOLUME, gamma);

    switch(interfaceAdvectionMethod_)
    {
    case CICSAM:

        //gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) + ib_.eqns(gamma) == 0.);
        gammaEqn_ = (cicsam::div(u, gradGamma, ib_.ibObjs(), gamma, timeStep, cicsam::HC) == 0.);
        break;
    case PLIC:

        gammaEqn_ = (plic::div(u, gradGamma, gamma, timeStep, geometries()["plicPolygons"]) == 0.);
        break;
    }

    return gammaEqn_.solve();
}

void PisoMultiphase::rhieChowInterpolation()
{
    Piso::rhieChowInterpolation();
    const Scalar sigma = surfaceTensionForce_->sigma();

    for(const Face& face: u.grid.interiorFaces())
    {
        const size_t id = face.id();
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar df = d.faces()[id];
        const Scalar rhof = rho.faces()[id];
        const Scalar rhoP = rho(lCell);
        const Scalar rhoQ = rho(rCell);

        Scalar kf;

        if(!u.grid.fluidCells().isInGroup(lCell))
            kf = surfaceTensionForce_->kappa()(rCell);
        else if(!u.grid.fluidCells().isInGroup(rCell))
            kf = surfaceTensionForce_->kappa()(lCell);
        else
            kf = surfaceTensionForce_->kappa().faces()[id];

        const Scalar gP = gamma(lCell);
        const Scalar gQ = gamma(rCell);

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        if(ib_.isIbCell(lCell) || ib_.isIbCell(rCell))
        {
            continue;
        }
        else
        {
            u(face) += sigma*df*kf*(gQ - gP)*rc/dot(rc, rc) - rhof*(g*d(lCell)*ft(lCell)/rhoP + (1. - g)*d(rCell)*ft(rCell)/rhoQ)/2.;
        }
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell& cellP = face.lCell();
        const Vector2D rf = face.centroid() - cellP.centroid();
        const Scalar df = d(face);
        const Scalar kf = surfaceTensionForce_->kappa()(face);
        const Scalar gf = gamma(face);
        const Scalar gP = gamma(cellP);
        const Scalar rhoP = rho(cellP);
        const Scalar rhof = rho(face);

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) += surfaceTensionForce_->sigma()*df*kf*(gf - gP)*rf/dot(rf, rf)
                    - rhof*d(cellP)*ft(cellP)/rhoP;
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            break;

        case VectorFiniteVolumeField::OUTFLOW:
            u(face) += dot(u(face), face.outwardNorm(cellP.centroid())) > 0. ? surfaceTensionForce_->sigma()*df*kf*(gf - gP)*rf/dot(rf, rf)
                                                                                                         - rhof*d(cellP)*ft(cellP)/rhoP : Vector2D(0., 0.);
            break;

        default:
            throw Exception("PisoMultiphase", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

//- External functions
