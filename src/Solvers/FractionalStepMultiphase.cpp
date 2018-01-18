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
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    //ib_.computeForce(rho, mu, u, p, g_);
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
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.5)
                 == ft.contactLineBcs(ib_));

    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update the gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    //- Update all other properties
    computeMomentumFlux(beta, timeStep);
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);

    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_.velocityBcs(u)
             == fv::laplacian(mu, u, 0.5) + src::src(ft0 + sg0, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar g = face.volumeWeight();
        u(face) = g * (u(face.lCell()) - timeStep / rho(face.lCell()) * (ft0(face.lCell()) + sg0(face.lCell())))
                  + (1. - g) * (u(face.rCell()) - timeStep / rho(face.rCell()) * (ft0(face.rCell()) + sg0(face.rCell())))
                  + timeStep / rho(face) * (ft(face) + sg(face));
    }

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &f: patch)
                {
                    const Cell &l = f.lCell();
                    u(f) = u(l) - timeStep / rho(l) * (ft0(l) + sg0(l))
                           + timeStep / rho(f) * (ft(f) + sg(f));
                }
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p, fluid_) + ib_.bcs(p) == src::div(u, fluid_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.computeFaces();
    gradP.faceToCell(rho, rho, fluid_);

    return error;
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho(cell) * (gradP(cell) - ft(cell) - sg(cell) + ft0(cell) + sg0(cell));

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) -= timeStep / rho(face) * gradP(face);
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }
}

void FractionalStepMultiphase::computeMomentumFlux(const ScalarFiniteVolumeField &beta, Scalar timeStep)
{
    rhoU.savePreviousTimeStep(timeStep, 1);

    rhoU.oldField(0).computeInteriorFaces([this, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta(f)) * gamma.oldField(0)(d) + beta(f) * gamma.oldField(0)(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.oldField(0).computeBoundaryFaces([this, &beta](const Face &f) {
        Scalar g = gamma.oldField(0)(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeInteriorFaces([this, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta(f)) * gamma(d) + beta(f) * gamma(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeBoundaryFaces([this](const Face &f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });
}

void FractionalStepMultiphase::updateProperties(Scalar timeStep)
{
    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell &cell) {
        Scalar g = gamma(cell);
        return (1. - g) * rho1_ + g * rho2_;
    });

    rho.computeFaces([this](const Face &face) {
        Scalar g = gamma(face);
        return (1. - g) * rho1_ + g * rho2_;
    });

    //grid_->sendMessages(rho); //- For correct gradient computation

    //- Update the gravitational source term
    gradRho.computeFaces();
    sg.savePreviousTimeStep(timeStep, 1.);
    for (const Face &face: grid().faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);

    sg.oldField(0).faceToCell(rho, rho.oldField(0), fluid_);
    sg.faceToCell(rho, rho, fluid_);

    //- Update viscosity from kinematic viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell &cell) {
        Scalar g = gamma(cell);
        return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });
    mu.computeFaces([this](const Face &face) {
        Scalar g = gamma(face);
        return rho(face) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    //- Update the surface tension
    ft.savePreviousTimeStep(timeStep, 1);
    ft.compute(ib_);

    //- Predicate ensures cell-centred values aren't overwritten for cells neighbouring ib cells
    auto p = [this](const Cell& cell)
    {
        for(const CellLink& nb: cell.neighbours())
            if(ib_.ibObj(nb.cell().centroid()))
                return false;
        return true;
    };

    ft.oldField(0).faceToCell(rho, rho.oldField(0), fluid_, p);
    ft.faceToCell(rho, rho, fluid_, p);
}
