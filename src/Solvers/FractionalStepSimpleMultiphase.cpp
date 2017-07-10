#include "FractionalStepSimpleMultiphase.h"
#include "ScalarGradient.h"
#include "Source.h"
#include "CrankNicolson.h"
#include "Cicsam.h"
#include "SeoMittal.h"

FractionalStepSimpleMultiphase::FractionalStepSimpleMultiphase(const Input &input,
                                                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStepSimple(input, grid),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gamma(addScalarField(input, "gamma")),
        rhoU(addVectorField("rhoU")),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gammaEqn_(input, gamma, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);
    g_ = input.caseInput().get<std::string>("Properties.g", "(0,0)");

    rho.copyBoundaryTypes(gamma);
    mu.copyBoundaryTypes(gamma);

    ft_ = std::make_shared<Celeste>(input, ib_, gamma, rho, mu, u, gradGamma);
    registerField(ft_);

    Scalar sigma = ft_->sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * sigma)));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);
}

void FractionalStepSimpleMultiphase::initialize()
{
    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);

    //- Set the fixed boundaries for rho and mu
    rho.setBoundaryFaces(ScalarFiniteVolumeField::FIXED, [this](const Face& face){
        Scalar g = gamma(face);
        return (1. - g)*rho1_ + g*rho2_;
    });

    mu.setBoundaryFaces(ScalarFiniteVolumeField::FIXED, [this](const Face& face){
        Scalar g = gamma(face);
        return (1. - g)*mu1_ + g*mu2_;
    });

    gamma.savePreviousTimeStep(0, 1);

    updateProperties(0.);
}

Scalar FractionalStepSimpleMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min(FractionalStepSimple::computeMaxTimeStep(maxCo, prevTimeStep), capillaryTimeStep_);
}

Scalar FractionalStepSimpleMultiphase::solve(Scalar timeStep)
{
    printf("Solving VoF equation...\n");
    solveGammaEqn(timeStep);

    printf("Solving momentum equation...\n");
    solveUEqn(timeStep);

    ib_.clearFreshCells();

    printf("Solving pressure eqn...\n");
    solvePEqn(timeStep);

    printf("Performing projection...\n");
    correctVelocity(timeStep);

    ib_.computeForce(rho, mu, u, p);

    ib_.update(timeStep);

    printf("Max divergence error = %.4e\n", grid_->comm().max(seo::maxDivergence(ib_, u)));
    printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

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

    grid_->sendMessages(gamma);

    gamma.interpolateFaces([this, &beta](const Face& face){
        return dot(u(face), face.outwardNorm(face.lCell().centroid())) > 0 ? 1. - beta(face): beta(face);
    });

    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepSimpleMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.bcs(u) == cn::laplacian(mu, u, 0.5)
                                                                          + fv::source(*ft_, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    u.interpolateFaces([](const Face &f){
        return f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
    });

    return error;
}

Scalar FractionalStepSimpleMultiphase::solvePEqn(Scalar timeStep)
{
    //pEqn_ = (seo::laplacian(ib_, rho, timeStep, p) == seo::div(ib_, u));
    pEqn_ = (fv::laplacian(timeStep/rho, p) + ib_.bcs(p) == fv::src::div(u));

    Scalar error = pEqn_.solve();

    grid_->sendMessages(p);

    p.interpolateFaces([](const Face& f){
        return f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
    });

    return error;
}

void FractionalStepSimpleMultiphase::correctVelocity(Scalar timeStep)
{
    //seo::correct(ib_, rho, p, *ft_, gradP, u, timeStep);
    gradP.computeWeighted(rho, fluid_);

    for(const Face& face: grid_->interiorFaces())
        u(face) -= timeStep/rho(face) * gradP(face);

    u.setBoundaryFaces(VectorFiniteVolumeField::NORMAL_GRADIENT, [this, timeStep](const Face& face) {
        return u(face) - timeStep/rho(face) * gradP(face);
    });

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep / rho(cell) * gradP(cell);
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
    rho.interpolateFaces(alpha);

    //- Update viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.compute([this](const Cell& cell){
        Scalar g = std::max(std::min(gamma(cell), 1.), 0.);
        return (1. - g)*mu1_ + g*mu2_;
    });
    mu.interpolateFaces(alpha);

    //- Update surface tension force
    ft_->savePreviousTimeStep(timeStep, 1);
    ft_->constructMatrices();
    ft_->compute();
}
