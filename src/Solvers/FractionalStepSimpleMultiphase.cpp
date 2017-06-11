#include <Discretization/GradientEvaluation.h>
#include "FractionalStepSimpleMultiphase.h"
#include "SeoMittal.h"
#include "Cicsam.h"
#include "Source.h"

FractionalStepSimpleMultiphase::FractionalStepSimpleMultiphase(const Input &input,
                                                               const Communicator &comm,
                                                               FiniteVolumeGrid2D &grid)
        :
        FractionalStepSimple(input,
                             comm,
                             grid),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gamma(addScalarField(input, "gamma")),
        rhoU(addVectorField("rhoU")),
        gradGamma(addVectorField("gradGamma")),
        gammaEqn_(input, comm, gamma, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);

    rho.copyBoundaryTypes(gamma);
    mu.copyBoundaryTypes(gamma);
}

void FractionalStepSimpleMultiphase::initialize()
{
    rho.compute([this](const Cell &c){
                    Scalar g = std::max(std::min(gamma(c), 1.), 0.);
                    return (1. - g)*rho1_ + g*rho2_;
                },
                [this](const Face &f){
                    Scalar g = std::max(std::min(gamma(f), 1.), 0.);
                    return (1. - g)*rho1_ + g*rho2_;
                });

    mu.compute([this](const Cell &c){
                    Scalar g = std::max(std::min(gamma(c), 1.), 0.);
                    return (1. - g)*mu1_ + g*mu2_;
                },
                [this](const Face &f){
                    Scalar g = std::max(std::min(gamma(f), 1.), 0.);
                    return (1. - g)*mu1_ + g*mu2_;
                });


    rho.savePreviousTimeStep(0, 1);
    mu.savePreviousTimeStep(0, 1);
}

Scalar FractionalStepSimpleMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);

    ib_.clearFreshCells();

    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_.computeForce(ScalarFiniteVolumeField(grid_, "rho", rho_), ScalarFiniteVolumeField(grid_, "mu", mu_), u, p);
    ib_.update(timeStep);

    comm_.printf("Max divergence error = %.4e\n", comm_.max(seo::maxDivergence(ib_, u)));
    comm_.printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

//- Private methods

Scalar FractionalStepSimpleMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    cicsam::interpolateFaces(u, gradGamma, gamma, timeStep, 0.5);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, gamma) +
                 ib_.bcs(gamma) == 0.);

    Scalar error = gammaEqn_.solve();

    rhoU.compute([this](const Face& f){
        Scalar g = std::max(std::min(gamma(f), 1.), 0.);
        return ((1. - g)*rho1_ + g*rho2_)*u(f);
    });

    rho.compute([this](const Cell &c){
        Scalar g = std::max(std::min(gamma(c), 1.), 0.);
        return (1. - g)*rho1_ + g*rho2_;
    });

    mu.compute([this](const Cell &c){
        Scalar g = std::max(std::min(gamma(c), 1.), 0.);
        return (1. - g)*mu1_ + g*mu2_;
    });

    auto alpha = [](const Face &f){
        Scalar l1 = (f.centroid() - f.lCell().centroid()).mag();
        Scalar l2 = (f.centroid() - f.rCell().centroid()).mag();
        return l2/(l1 + l2);
    };

    rho.interpolateFaces(alpha);
    mu.interpolateFaces(alpha);

    return error;
}

Scalar FractionalStepSimpleMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.bcs(u) == fv::laplacian(mu, u));
    Scalar error = uEqn_.solve();

    u.interpolateFaces([](const Face &f){
        Scalar l1 = (f.centroid() - f.lCell().centroid()).mag();
        Scalar l2 = (f.centroid() - f.rCell().centroid()).mag();
        return l2/(l1 + l2);
    });

    return error;
}

Scalar FractionalStepSimpleMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep/rho, p) + ib_.bcs(p) == source::div(u));
    Scalar error = pEqn_.solve();

    p.interpolateFaces([](const Face& f){
        Scalar l1 = (f.lCell().centroid() - f.centroid()).mag();
        Scalar l2 = (f.rCell().centroid() - f.centroid()).mag();

        return l2/(l1 + l2);
    });

    return error;
}

void FractionalStepSimpleMultiphase::correctVelocity(Scalar timeStep)
{
    fv::computeGradient(fv::FACE_TO_CELL, fluid_, p, gradP);
    for(const Face& face: grid_.interiorFaces())
        u(face) -= timeStep/rho(face)*gradP(face);

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep/rho(cell)*gradP(cell);

    u.setBoundaryFaces();

    //seo::correct(ib_, rho, p, gradP, u, timeStep);
}
