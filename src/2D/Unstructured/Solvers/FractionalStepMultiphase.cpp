#include "Math/Algorithm.h"
#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/SecondOrderExplicitDivergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/Discretization/Cicsam.h"

#include "FractionalStepMultiphase.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStep(input, grid),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      gamma_(*addField<Scalar>(input, "gamma", fluid_)),
      beta_(*addField<Scalar>("beta", fluid_)),
      rhoU_(*addField<Vector2D>("rhoU", fluid_)),
      fst_(input, grid_, fluid_),
      sg_(*addField<Vector2D>("sg", fluid_)),
      gradGamma_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(gamma_, fluid_)))),
      gradRho_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(rho_, fluid_)))),
      gammaEqn_(input, gamma_, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", FractionalStep::rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", FractionalStep::rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", FractionalStep::mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", FractionalStep::mu_);

    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();
    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4. * M_PI * fst_.sigma())));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);

    addField(fst_.fst());
    addField(fst_.gammaTilde());
    addField(fst_.kappa());
    addField<Vector2D>(fst_.gradGammaTilde());
    addField(fst_.n());
}

void FractionalStepMultiphase::initialize()
{
    FractionalStep::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma_.compute(*fluid_);
    u_.savePreviousTimeStep(0, 2);
    gamma_.savePreviousTimeStep(0, 2);
    gradGamma_.savePreviousTimeStep(0, 2);
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

    printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

//- Private methods

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    //- Advect volume fractions
    gamma_.savePreviousTimeStep(timeStep, 2);
    gammaEqn_ = (fv::ddt(gamma_, timeStep) + cicsam::div(u_, gamma_, gradGamma_, timeStep, 0.5)
                 == 0.);

    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma_);
    gamma_.interpolateFaces();

    //- Momentum flux
    cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_, gradGamma_, timeStep, rhoU_);

    //- Update the gradient
    gradGamma_.savePreviousTimeStep(timeStep, 2);
    gradGamma_.compute(*fluid_);
    grid_->sendMessages(gradGamma_);

    //- Update all other properties
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 2);
    const auto &fst = *fst_.fst();

    uEqn_ = (fv::ddt(rho_, u_, timeStep) + fv::div2e(rhoU_, u_, 0.5)
             == fv::laplacian(mu_, u_, 0.5) + src::src(fst + sg_));

    Scalar error = uEqn_.solve();

    grid_->sendMessages(u_);

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar g = face.volumeWeight();
        u_(face) = g * (u_(face.lCell()) - timeStep / rho_(face.lCell()) * (fst(face.lCell()) + sg_(face.lCell())))
                + (1. - g) * (u_(face.rCell()) - timeStep / rho_(face.rCell()) * (fst(face.rCell()) + sg_(face.rCell())))
                + timeStep / rho_(face) * (fst(face) + sg_(face));
    }

    for (const FaceGroup &patch: grid_->patches())
        switch (u_.boundaryType(patch))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            for (const Face &f: patch)
            {
                const Cell &l = f.lCell();
                u_(f) = u_(l) - timeStep / rho_(l) * (fst(l) + sg_(l))
                        + timeStep / rho_(f) * (fst(f) + sg_(f));
            }
            break;
        case VectorFiniteVolumeField::SYMMETRY:
            for (const Face &f: patch)
                u_(f) = u_(f.lCell()) - dot(u_(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
            break;
        }

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p_) == src::div(u_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p_);

    p_.setBoundaryFaces();
    gradP_.computeFaces();
    gradP_.faceToCell(rho_, rho_, *fluid_);

    return error;
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: *fluid_)
        u_(cell) -= timeStep / rho_(cell) * gradP_(cell);

    grid_->sendMessages(u_);

    for (const Face &face: grid_->interiorFaces())
        u_(face) -= timeStep / rho_(face) * gradP_(face);

    for (const FaceGroup &patch: grid_->patches())
        switch (u_.boundaryType(patch))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            for (const Face &face: patch)
                u_(face) -= timeStep / rho_(face) * gradP_(face);
            break;
        case VectorFiniteVolumeField::SYMMETRY:
            for (const Face &f: patch)
                u_(f) = u_(f.lCell()) - dot(u_(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
            break;
        }
}

void FractionalStepMultiphase::updateProperties(Scalar timeStep)
{
    //- Update density
    rho_.savePreviousTimeStep(timeStep, 1);

    rho_.computeCells([this](const Cell &cell) {
        Scalar g = gamma_(cell);
        return (1. - g) * rho1_ + g * rho2_;
    });

    rho_.computeFaces([this](const Face &face) {
        Scalar g = gamma_(face);
        return (1. - g) * rho1_ + g * rho2_;
    });

    //- Update the gravitational source term
    gradRho_.computeFaces();

    for (const Face &face: grid_->faces())
        sg_(face) = dot(g_, -face.centroid()) * gradRho_(face);

    sg_.faceToCell(rho_, rho_, *fluid_);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(sg_);

    //- Update viscosity from kinematic viscosity
    mu_.savePreviousTimeStep(timeStep, 1);

    mu_.computeCells([this](const Cell &cell) {
        Scalar g = gamma_(cell);
        return (1. - g) * mu1_ + g * mu2_;
        //return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    mu_.computeFaces([this](const Face &face) {
        Scalar g = gamma_(face);
        return (1. - g) * mu1_ + g * mu2_;
        //return rho(face) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    //- Update the surface tension
    fst_.computeFaceInterfaceForces(gamma_, gradGamma_);
    fst_.fst()->faceToCell(rho_, rho_, *fluid_);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(*fst_.fst());
}
