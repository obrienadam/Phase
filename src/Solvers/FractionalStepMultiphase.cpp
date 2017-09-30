#include "FractionalStepMultiphase.h"
#include "Source.h"
#include "Cicsam.h"
#include "Hric.h"
#include "SeoMittal.h"
#include "Algorithm.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input,
                                                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gamma(addScalarField(input, "gamma")),
        rhoU(addVectorField("rhoU")),
        ft(addVectorField(std::make_shared<Celeste>(input, ib_, gamma, rho, mu, u, gradGamma))),
        sg(addVectorField("sg")),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gradRho(addVectorField(std::make_shared<ScalarGradient>(rho))),
        gammaEqn_(input, gamma, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);

    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();
    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * ft.sigma())));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);
}

void FractionalStepMultiphase::initialize()
{
    FractionalStep::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);
    u.savePreviousTimeStep(0, 1);
    gamma.savePreviousTimeStep(0, 1);
    updateProperties(0.);
    updateProperties(0.);
}

Scalar FractionalStepMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min(FractionalStep::computeMaxTimeStep(maxCo, prevTimeStep), capillaryTimeStep_);
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);

    ib_.clearFreshCells();

    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_.computeForce(rho, mu, u, p);
    ib_.update(timeStep);

    printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

//- Private methods

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    auto beta = cicsam::beta(u, gradGamma, gamma, timeStep, 0.5);

    //- Advect volume fractions
    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.5) == 0.);

    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma);

    gamma.interpolateFaces();

    //- Update mass fluxes
    rhoU.savePreviousTimeStep(timeStep, 1);

    rhoU.oldField(0).computeInteriorFaces([this, &beta](const Face& f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell& d = flux > 0. ? f.lCell(): f.rCell();
        const Cell& a = flux > 0. ? f.rCell(): f.lCell();
        Scalar g = (1. - beta(f)) * gamma.oldField(0)(d) + beta(f) * gamma.oldField(0)(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.oldField(0).computeBoundaryFaces([this, &beta](const Face& f) {
        Scalar g = gamma.oldField(0)(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeInteriorFaces([this, &beta](const Face& f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell& d = flux > 0. ? f.lCell(): f.rCell();
        const Cell& a = flux > 0. ? f.rCell(): f.lCell();
        Scalar g = (1. - beta(f)) * gamma(d) + beta(f) * gamma(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeBoundaryFaces([this](const Face& f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    //- Update the gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    //- Update all other properties
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_.bcs(u)
             == fv::laplacian(mu, u, 0.5));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    u.interpolateFaces();

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib_.bcs(p) == src::div(u));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.computeFaces();

    return error;
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    u.setBoundaryFaces(VectorFiniteVolumeField::NORMAL_GRADIENT, [this, timeStep](const Face &face) {
        return u(face) - timeStep / rho(face) * gradP(face);
    });

    u.setBoundaryFaces(VectorFiniteVolumeField::SYMMETRY, [this, timeStep](const Face &face) {
        const Vector2D &nw = face.norm();
        return u(face.lCell()) - dot(u(face.lCell()), nw) * nw / nw.magSqr();
    });

    gradP.faceToCell(rho, rho, fluid_);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho(cell) * gradP(cell);
}

void FractionalStepMultiphase::updateProperties(Scalar timeStep)
{
    //- Update density, density gradient and gravity source
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell& c) {
        Scalar g = gamma(c);
        return (1. - g) * rho1_ + g * rho2_;
    });

    grid_->sendMessages(rho);

    rho.computeFaces([this](const Face& f) {
        Scalar g = gamma(f);
        return (1. - g) * rho1_ + g * rho2_;
    });

    gradRho.computeFaces();

    sg.fill(Vector2D(0., 0.));
    for (const Face &face: grid_->faces())
        sg(face) = -dot(g_, face.centroid()) * gradRho(face);

    sg.faceToCell(rho, rho, fluid_);

    //- Update viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell &cell) {
        Scalar g = clamp(gamma(cell), 0., 1.);
        return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    grid_->sendMessages(mu);

    mu.computeFaces([this](const Face& f) {
        Scalar g = gamma(f);
        return (1. - g) * mu1_ + g * mu2_;
    });

    mu.interpolateFaces();

    //- Update surface tension force
    ft.savePreviousTimeStep(timeStep, 1);
    //ft_->constructMatrices();
    //ft_->compute();
}
