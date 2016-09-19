#include "PisoMultiphase.h"
#include "Cicsam.h"
#include "Plic.h"
#include "Celeste.h"
#include "CrankNicolson.h"
#include "FaceInterpolation.h"

PisoMultiphase::PisoMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
      ft(addVectorField("ft")),
      gammaEqn_(input, gamma, "gamma")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    //- Configuration
    interfaceAdvectionMethod_ = CICSAM;
    const std::string tmp = input.caseInput().get<std::string>("Solver.surfaceTensionModel");

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

    Scalar sigma = surfaceTensionForce_->sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for(const Face& face: grid_.interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_, sqrt(((rho1_ + rho2_)*delta*delta*delta)/(4*M_PI*sigma)));
    }

    printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

Scalar PisoMultiphase::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for(size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();
        solveUEqn(timeStep);

        for(size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            solvePCorrEqn();
            correctPressure();
            correctVelocity();
        }
    }

    solveGammaEqn(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.; // just to get rid of warning
}

Scalar PisoMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min(
                Piso::computeMaxTimeStep(maxCo, prevTimeStep),
                capillaryTimeStep_
                );
}

//- Protected methods

void PisoMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &alpha = surfaceTensionForce_->gammaTilde();

    rho.savePreviousTimeStep(0, 1);

    for(const Cell& cell: rho.grid.activeCells())
        rho(cell) = (1. - alpha(cell))*rho1_ + alpha(cell)*rho2_;

    harmonicInterpolateFaces(fv::INVERSE_VOLUME, rho);
}

void PisoMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &alpha = surfaceTensionForce_->gammaTilde();

    mu.savePreviousTimeStep(0, 1);

    for(const Cell& cell: mu.grid.activeCells())
        mu(cell) = (1. - alpha(cell))*mu1_ + alpha(cell)*mu2_;

    interpolateFaces(fv::INVERSE_VOLUME, mu);
}

Scalar PisoMultiphase::solveUEqn(Scalar timeStep)
{
    ft = surfaceTensionForce_->compute(); // surface tension force. The order matters here since the smoothed gamma field is used for rho and mu
    computeRho();
    computeMu();

    sg = fv::gravity(rho, g_);

    for(const Cell& cell: gamma.grid.fluidCells())
        ft(cell) *= 2.*rho(cell)/(rho1_ + rho2_);

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
        gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::cn(u, gradGamma, surfaceTensionForce_->n(), gamma, timeStep) + ib_.eqns(gamma) == 0.);
        //gammaEqn_ = (fv::ddt(gamma, timeStep) + fv::div(u, gamma) + ib_.eqns(gamma) == 0.);
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

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Scalar df = d(face);
        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        if(!(grid_.fluidCells().isInGroup(lCell) && grid_.fluidCells().isInGroup(rCell)))
            continue;

        u(face) += df*ft(face) - (g*d(lCell)*ft(lCell) + (1. - g)*d(rCell)*ft(rCell));
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell& cellP = face.lCell();
        const Scalar df = d(face);

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) += df*ft(face) - d(cellP)*ft(cellP);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            break;

        default:
            throw Exception("PisoMultiphase", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

//- External functions
