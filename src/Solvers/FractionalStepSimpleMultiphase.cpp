#include "FractionalStepSimpleMultiphase.h"
#include "Source.h"
#include "CrankNicolson.h"
#include "Cicsam.h"
#include "Hric.h"
#include "SeoMittal.h"
#include "Algorithm.h"

FractionalStepSimpleMultiphase::FractionalStepSimpleMultiphase(const Input &input,
                                                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStepSimple(input, grid),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gamma(addScalarField(input, "gamma")),
        rhoU(addVectorField("rhoU")),
        ft(addVectorField(std::make_shared<Celeste>(input, ib_, gamma, rho, mu, u, gradGamma))),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gammaEqn_(input, gamma, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);
    g_ = input.caseInput().get<std::string>("Properties.g", "(0,0)");

    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();
    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * ft.sigma())));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);
}

void FractionalStepSimpleMultiphase::initialize()
{
    FractionalStepSimple::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);
    gamma.savePreviousTimeStep(0, 1);
    updateProperties(0.);
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
    auto beta = cicsam::beta(u, gradGamma, gamma, timeStep, 0.5);

    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.) + ib_.bcs(gamma) == 0.);
    Scalar error = gammaEqn_.solve();

    //- Send messages and update faces
    grid_->sendMessages(gamma);

    //- Update gamma face
    gamma.oldField(0).computeInteriorFaces([this, &beta](const Face& f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell& d = flux > 0. ? f.lCell(): f.rCell();
        const Cell& a = flux > 0. ? f.rCell(): f.lCell();
        return (1. - beta(f)) * gamma.oldField(0)(d) + beta(f) * gamma.oldField(0)(a);
    });
    gamma.oldField(0).setBoundaryFaces();

    gamma.computeInteriorFaces([this, &beta](const Face& f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell& d = flux > 0. ? f.lCell(): f.rCell();
        const Cell& a = flux > 0. ? f.rCell(): f.lCell();
        return (1. - beta(f)) * gamma(d) + beta(f) * gamma(a);
    });
    gamma.setBoundaryFaces();

    //- Update the gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepSimpleMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.) + ib_.bcs(u)
             == fv::laplacian(mu, u, 0.5) + src::src(ft, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    u.interpolateFaces([](const Face &f){
        return f.volumeWeight();
    });

    return error;
}

Scalar FractionalStepSimpleMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep/rho, p) + ib_.bcs(p) == src::div(u));
    Scalar error = pEqn_.solve();

    grid_->sendMessages(p);
    p.setBoundaryFaces();
    gradP.compute(rho);

    return error;
}

void FractionalStepSimpleMultiphase::correctVelocity(Scalar timeStep)
{
    for(const Face& face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    u.setBoundaryFaces(VectorFiniteVolumeField::NORMAL_GRADIENT, [this, timeStep](const Face& face) {
        return u(face) - timeStep/rho(face) * gradP(face);
    });

    u.setBoundaryFaces(VectorFiniteVolumeField::SYMMETRY, [this, timeStep] (const Face& face) {
        const Vector2D &nw = face.norm();
        return u(face.lCell()) - dot(u(face.lCell()), nw) * nw / nw.magSqr();
    });

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep / rho(cell) * gradP(cell);
}

void FractionalStepSimpleMultiphase::updateProperties(Scalar timeStep)
{
    //- Update mass flux, make sure time steps are synced with gamma
    rhoU.computeFaces([this](const Face& f) {
        Scalar g = gamma.oldField(0)(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });
    rhoU.savePreviousTimeStep(timeStep, 1);
    rhoU.computeFaces([this](const Face& f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell& cell){
        Scalar g = clamp(gamma(cell), 0., 1.);
        return (1. - g)*rho1_ + g*rho2_;
    });
    rho.interpolateFaces([](const Face& f) { return f.volumeWeight(); });
    rho.computeBoundaryFaces([this](const Face& f) {
        Scalar g = clamp(gamma(f), 0., 1.);
        return (1. - g)*rho1_ + g*rho2_;
    });

    //- Update viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell& cell){
        Scalar g = clamp(gamma(cell), 0., 1.);
        return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });
    mu.computeFaces([this](const Face& f) {
        Scalar g = clamp(gamma(f), 0., 1.);
        return rho(f) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    //- Update surface tension force
    ft.savePreviousTimeStep(timeStep, 1);
    //ft_->constructMatrices();
    //ft_->compute();
}
