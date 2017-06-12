#include <Discretization/GradientEvaluation.h>
#include "FractionalStepSimpleMultiphase.h"
#include "SeoMittal.h"
#include "Cicsam.h"
#include "Source.h"
#include "CrankNicolson.h"

FractionalStepSimpleMultiphase::FractionalStepSimpleMultiphase(const Input &input,
                                                               const Communicator &comm,
                                                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
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
    //- Ensure the computation starts with a valid gamma field
    fv::computeGradient(fv::FACE_TO_CELL, fluid_, gamma, gradGamma);

    //- Be careful about this, mu and rho must have valid time histories!!!
    rho.computeBoundaryFaces([this](const Face& face){
        Scalar g = gamma(face);
        return (1. - g)*rho1_ + g*rho2_;
    });

    mu.computeBoundaryFaces([this](const Face& face){
        Scalar g = gamma(face);
        return (1. - g)*mu1_ + g*mu2_;
    });

    gamma.savePreviousTimeStep(0, 1);

    updateProperties(0.);
}

Scalar FractionalStepSimpleMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);

    ib_.clearFreshCells();

    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_.computeForce(ScalarFiniteVolumeField(grid(), "rho", rho_),
                     ScalarFiniteVolumeField(grid(), "mu", mu_), u, p);

    ib_.update(timeStep);

    comm_.printf("Max divergence error = %.4e\n", comm_.max(seo::maxDivergence(ib_, u)));
    comm_.printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

//- Private methods

Scalar FractionalStepSimpleMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    auto beta = cicsam::computeBeta(u, gradGamma, gamma, timeStep, 0.5);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma) +
                 ib_.bcs(gamma) == 0.);

    Scalar error = gammaEqn_.solve();

    // - While this may affect mass conservation, it prevents issues at high density ratios
    gamma.compute([this](const Cell &cell){
        return std::max(std::min(gamma(cell), 1.), 0.);
    });

    grid_->sendMessages(comm_, gamma);

    gamma.interpolateFaces([this, &beta](const Face& face){
        return dot(u(face), face.outwardNorm(face.lCell().centroid())) > 0 ? 1. - beta(face): beta(face);
    });

    gamma.setBoundaryFaces();
    fv::computeGradient(fv::FACE_TO_CELL, fluid_, gamma, gradGamma);

    grid_->sendMessages(comm_, gradGamma);

    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepSimpleMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.bcs(u) == cn::laplacian(mu, u, 0.5));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(comm_, u);

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

    grid_->sendMessages(comm_, p);

    p.interpolateFaces([](const Face& f){
        Scalar l1 = (f.lCell().centroid() - f.centroid()).mag();
        Scalar l2 = (f.rCell().centroid() - f.centroid()).mag();

        return l2/(l1 + l2);
    });

    fv::computeInverseWeightedGradient(rho, p, gradP);

    grid_->sendMessages(comm_, gradP);

    return error;
}

void FractionalStepSimpleMultiphase::correctVelocity(Scalar timeStep)
{
    for(const Face& face: grid_->interiorFaces())
        u(face) -= timeStep/rho(face)*gradP(face);

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep/rho(cell)*gradP(cell);

    for (const Patch& patch: grid_->patches())
    {
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face &face: patch)
                {
                    Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
                }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face& face: patch)
                    u(face) -= timeStep / rho(face) * gradP(face);
                break;
        };
    }

    //seo::correct(ib_, rho, p, gradP, u, timeStep);
}

void FractionalStepSimpleMultiphase::updateProperties(Scalar timeStep)
{
    auto alpha = [](const Face& f){
        return f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
    };

    //- Update rhoU
    rhoU.compute([this](const Face& face){
        Scalar g = std::max(std::min(gamma(face), 1.), 0.);
        Scalar g0 = std::max(std::min(gamma.oldField(0)(face), 1.), 0.);
        Scalar rhoF = (1. - g)*rho1_ + g*rho2_;
        Scalar rhoF0 = (1. - g0)*rho1_ + g0*rho2_;
        return (rhoF + rhoF0)/2.*u(face);
    });

    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.compute([this](const Cell& cell){
        Scalar g = std::max(std::min(gamma(cell), 1.), 0.);
        return (1. - g)*rho1_ + g*rho2_;
    });
    grid_->sendMessages(comm_, rho);
    rho.interpolateFaces(alpha);

    //- Update viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.compute([this](const Cell& cell){
        Scalar g = std::max(std::min(gamma(cell), 1.), 0.);
        return (1. - g)*mu1_ + g*mu2_;
    });
    grid_->sendMessages(comm_, mu);
    mu.interpolateFaces(alpha);
}