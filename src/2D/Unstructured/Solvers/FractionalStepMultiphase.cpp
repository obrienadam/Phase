#include "Math/Algorithm.h"
#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/Discretization/Cicsam.h"
#include "FiniteVolume/Multiphase/Celeste.h"

#include "FractionalStepMultiphase.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid),
        rho(*addField<Scalar>("rho")),
        mu(*addField<Scalar>("mu")),
        gamma(*addField<Scalar>(input, "gamma")),
        beta(*addField<Scalar>("beta")),
        rhoU(*addField<Vector2D>("rhoU")),
        ft(*std::static_pointer_cast<Celeste>(addField<Vector2D>(std::make_shared<Celeste>(input, grid_, ib_)))),
        sg(*addField<Vector2D>("sg")),
        gradGamma(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(gamma)))),
        gradRho(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(rho)))),
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

    addField<Scalar>(ft.gammaTilde());
    addField<Scalar>(ft.kappa());
    addField<Vector2D>(ft.gradGammaTilde());
    addField<Vector2D>(ft.n());
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

    //ib_->computeForce(rho, mu, u, p, g_);
    ib_->update(timeStep);

    printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

//- Private methods

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    cicsam::beta(u, gradGamma, gamma, timeStep, beta, 0.5);

    //- Advect volume fractions
    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.5)
                 == ft.contactLineBcs(gamma));

    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update the gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    //- Update all other properties
    cicsam::computeMomentumFlux(rho1_, rho2_, u, gamma, beta, timeStep, rhoU);
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_->velocityBcs(u)
             == fv::laplacian(mu, u, 0.5) + src::src(ft + sg));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar g = face.volumeWeight();
        u(face) = g * (u(face.lCell()) - timeStep / rho(face.lCell()) * (ft(face.lCell()) + sg(face.lCell())))
                  + (1. - g) * (u(face.rCell()) - timeStep / rho(face.rCell()) * (ft(face.rCell()) + sg(face.rCell())))
                  + timeStep / rho(face) * (ft(face) + sg(face));
    }

    for (const FaceGroup &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &f: patch)
                {
                    const Cell &l = f.lCell();
                    u(f) = u(l) - timeStep / rho(l) * (ft(l) + sg(l))
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
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib_->bcs(p) == src::div(u));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.computeFaces();
    gradP.faceToCell(rho, rho, fluid_);

    return error;
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho(cell) * gradP(cell);

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    for (const FaceGroup &patch: grid_->patches())
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

    //- Update the gravitational source term
    gradRho.computeFaces();
    for (const Face &face: grid_->faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);

    sg.faceToCell(rho, rho, fluid_);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(sg);

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
    ft.computeFaceInterfaceForces(gamma, gradGamma);
    ft.faceToCell(rho, rho, fluid_, p);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(ft);
}
