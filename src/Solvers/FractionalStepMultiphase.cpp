#include "FractionalStepMultiphase.h"
#include "CrankNicolson.h"
#include "Cicsam.h"
#include "ScalarGradient.h"
#include "FaceInterpolation.h"
#include "Source.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input,
                                                   std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid),
        gamma(addScalarField(input, "gamma")),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gradRho(addVectorField(std::make_shared<ScalarGradient>(rho))),
        sg(addVectorField("sg")),
        rhoU(addVectorField("rhoU")),
        gammaEqn_(input, gamma, "gammaEqn")
{
    cicsamBlending_ = input.caseInput().get<Scalar>("Solver.cicsamBlending", 1.);
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);

    rho.copyBoundaryTypes(gamma);
    mu.copyBoundaryTypes(gamma);

    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    ft_ = std::make_shared<Celeste>(input, ib_, gamma, rho, mu, u, gradGamma);

    registerField(ft_);
    ft_->registerSubFields(*this);

    Scalar sigma = ft_->sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * sigma)));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);
    printf("CICSAM blending constant (k): %.2f\n", cicsamBlending_);
    printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

void FractionalStepMultiphase::initialize()
{
    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);

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

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep); //- Solve a sharp gamma equation
    solveUEqn(timeStep); //- Solve a momentum prediction

    solvePEqn(timeStep); //- Solve the pressure equation, using sharp value of rho
    correctVelocity(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    printf("Max absolute velocity divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

Scalar FractionalStepMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min( //- Both args have already been globalling minimized
            FractionalStep::computeMaxTimeStep(maxCo, prevTimeStep),
            capillaryTimeStep_
    );
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    auto beta = cicsam::computeBeta(u, gradGamma, gamma, timeStep, cicsamBlending_);
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

    //- Recompute gradGamma for next cicsam iteration
    gamma.setBoundaryFaces();
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma); //- In case donor cell is on another proc

    updateProperties(timeStep);

    return error;
}

//- Protected methods
Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.bcs(u)
             == cn::laplacian(mu, u, 0.5) - fv::source(gradP - sg.oldField(0) - ft_->oldField(0), fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);
    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, phi) + ib_.bcs(phi) == fv::src::div(u));

    Scalar error = pEqn_.solve();

    for(const Cell& cell: grid_->localActiveCells())
        p(cell) += phi(cell);

    grid_->sendMessages(p);

    p.interpolateFaces([](const Face& face){
        Scalar l1 = (face.centroid() - face.lCell().centroid()).mag();
        Scalar l2 = (face.centroid() - face.rCell().centroid()).mag();
        return l2/(l1 + l2);
    });

    //- Compute pressure gradient
    gradP.savePreviousTimeStep(timeStep, 1);

    //- Weighted gradients greatly reduce the effect of large pressure differences
    gradP.computeWeighted(rho, fluid_);
    grid_->sendMessages(gradP);

    return error;
}

void FractionalStepMultiphase::computeFaceVelocities(Scalar timeStep)
{
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft = *ft_;
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);

    for (const Face &face: u.grid().interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar l1 = (face.centroid() - face.lCell().centroid()).mag();
        Scalar l2 = (face.centroid() - face.rCell().centroid()).mag();
        Scalar g = l2/(l1 + l2);

        u(face) = g * (u(lCell) + timeStep / rho0(lCell) * (gradP(lCell) - sg0(lCell) - ft0(lCell)))
                  + (1. - g) * (u(rCell) + timeStep / rho0(rCell) * (gradP(rCell) - sg0(rCell) - ft0(rCell)))
                  + timeStep / rho(face) * (sg(face) + ft(face) - gradP(face));
    }

    for (const Face &face: u.grid().boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) = u(cell) + timeStep / rho0(cell) * (gradP(cell) - sg0(cell) - ft0(cell)) +
                          timeStep / rho(face) * (sg(face) + ft(face) - gradP(face));
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(cell.centroid());
                u(face) = u(cell) - dot(u(cell), nWall) * nWall / nWall.magSqr();
            }
                break;

                //default:
                //throw Exception("FractionalStep", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft = *ft_;
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep * ((gradP(cell) - sg(cell) - ft(cell)) / rho(cell) -
                               (gradP0(cell) - sg0(cell) - ft0(cell)) / rho0(cell));

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face));

    for (const Face &face: grid_->boundaryFaces())
    {
        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
            }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face));
                break;
        };
    }
}

void FractionalStepMultiphase::updateProperties(Scalar timeStep)
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

    //- Update the surface tension
    ft_->savePreviousTimeStep(timeStep, 1);
    ft_->compute();

    //- Update the gravitational source term
    gradRho.computeWeighted(rho, fluid_);

    sg.savePreviousTimeStep(timeStep, 1);
    sg.compute([this](const Cell& cell) {
        return dot(g_, -cell.centroid()) * gradRho(cell);
    }, [this](const Face& face) {
        return dot(g_, -face.centroid()) * gradRho(face);
    });
}

void FractionalStepMultiphase::checkMassFluxConsistency(Scalar timeStep)
{
    Scalar max = 0.;

    for (const Face &face: grid_->faces())
    {
        Scalar gammaMom = (rho(face) - rho1_) / (rho2_ - rho1_);
        Scalar gammaVof = gamma(face);

        Scalar error = fabs((gammaMom - gammaVof));

        if (error > max)
            max = error;
    }

    printf("Max mass flux error = %.2lf.\n", max);
}
