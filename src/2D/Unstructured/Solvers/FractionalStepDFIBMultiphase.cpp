#include "FractionalStepDFIBMultiphase.h"
#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/Discretization/Cicsam.h"

FractionalStepDirectForcingMultiphase::FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepDFIB(input, grid),
      gamma_(*addField<Scalar>(input, "gamma", fluid_)),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      sg_(*addField<Vector2D>("sg", fluid_)),
      rhoU_(grid_, "rhoU", Vector2D(0., 0.), true, false, fluid_),
      gradGamma_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(gamma_, fluid_)))),
      gradRho_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(rho_, fluid_)))),
      fst_(input, grid_, fluid_, ib_),
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
                                      std::sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * fst_.sigma())));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);

    addField(fst_.fst());
    addField(fst_.kappa());
    addField(fst_.gammaTilde());
    addField<Vector2D>(fst_.gradGammaTilde());
    addField(fst_.n());
}

void FractionalStepDirectForcingMultiphase::initialize()
{
    FractionalStepDFIB::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma_.compute(*fluid_);
    u_.savePreviousTimeStep(0, 1);

    updateProperties(0.);
    updateProperties(0.);
}

Scalar FractionalStepDirectForcingMultiphase::solve(Scalar timeStep)
{
    //- Perform field extension
    //solveExtEqns();

    grid_->comm().printf("Updating IB positions and cell categories...\n");
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    grid_->comm().printf("Solving gamma equation...\n");
    solveGammaEqn(timeStep);

    grid_->comm().printf("Updating physical properties...\n");
    updateProperties(timeStep);

    grid_->comm().printf("Solving momentum and pressure equations...\n");
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    grid_->comm().printf("Computing IB forces...\n");
    //ib_->computeForce(rho_, mu_, u_, p_, gamma_, fst_, g_);

    return 0;
}

Scalar FractionalStepDirectForcingMultiphase::solveGammaEqn(Scalar timeStep)
{
    auto beta = cicsam::beta(u_, gradGamma_, gamma_, timeStep, 0.5);

    //- Advect volume fractions
    gamma_.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma_, timeStep) + cicsam::div(u_, beta, gamma_, 0.) == 0.);
    fst_.contactLineBcs(gammaEqn_);

    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma_);
    gamma_.interpolateFaces();

    //- Update the gradient
    gradGamma_.compute(*fluid_);
    grid_->sendMessages(gradGamma_);

    //- Update all other properties
    cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_, beta, timeStep, rhoU_);

    return error;
}

Scalar FractionalStepDirectForcingMultiphase::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);
    const auto &fst = *fst_.fst();

    gradP_.faceToCell(rho_, rho_.oldField(0), *fluid_);

    uEqn_ = (fv::ddt(rho_, u_, timeStep) + fv::div(rhoU_, u_, 0.)
             == fv::laplacian(mu_, u_, 0.) + src::src(fst + sg_ - gradP_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u_);

    ib_->computeForcingTerm(rho_, u_, timeStep, fb_);

    for (const Cell &cell: grid_->cells())
        u_(cell) = u_.oldField(0)(cell);

    uEqn_ = (fv::ddt(rho_, u_, timeStep) + fv::div(rhoU_, u_, 0.)
             == fv::laplacian(mu_, u_, 0.) + src::src(fst + sg_ + fb_ - gradP_));

    error = uEqn_.solve();

    for(const Cell& cell: *fluid_)
        u_(cell) += timeStep / rho_(cell) * gradP_(cell);

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

Scalar FractionalStepDirectForcingMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p_) + ib_->bcs(p_) == src::div(u_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p_);

    p_.setBoundaryFaces();
    gradP_.computeFaces();
    gradP_.faceToCell(rho_, rho_, *fluid_);

    return error;
}

void FractionalStepDirectForcingMultiphase::updateProperties(Scalar timeStep)
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
    });

    mu_.computeFaces([this](const Face &face) {
        Scalar g = gamma_(face);
        return (1. - g) * mu1_ + g * mu2_;
    });

    //- Update the surface tension
    fst_.computeFaceInterfaceForces(gamma_, gradGamma_);
    fst_.fst()->faceToCell(rho_, rho_, *fluid_);
    fst_.fst()->fill(Vector2D(0., 0.), ib_->solidCells());

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(*fst_.fst());
}

void FractionalStepDirectForcingMultiphase::correctVelocity(Scalar timeStep)
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
