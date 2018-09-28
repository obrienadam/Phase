#include "Math/Algorithm.h"
#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricLaplacian.h"
#include "FiniteVolume/Discretization/AxisymmetricSource.h"
#include "FiniteVolume/Discretization/AxisymmetricStressTensor.h"
#include "FiniteVolume/Discretization/AxisymmetricCicsam.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"

#include "FractionalStepAxisymmetricDFIBMultiphase.h"

FractionalStepAxisymmetricDFIBMultiphase::FractionalStepAxisymmetricDFIBMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepAxisymmetricDFIB(input, grid),
      gamma_(*addField<Scalar>(input, "gamma", fluid_)),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      rhoU_(*addField<Vector2D>("rhoU", fluid_)),
      fst_(*addField<Vector2D>("fst", fluid_)),
      sg_(*addField<Vector2D>("sg", fluid_)),
      gradGamma_(static_cast<ScalarGradient&>(*addField<Vector2D>(std::make_shared<ScalarGradient>(gamma_, fluid_)))),
      gammaEqn_(input, gamma_, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", FractionalStep::rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", FractionalStep::rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", FractionalStep::mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", FractionalStep::mu_);
}

void FractionalStepAxisymmetricDFIBMultiphase::initialize()
{
    FractionalStepAxisymmetricDFIB::initialize();
    gradGamma_.compute(*fluid_);
    updateProperties(0.);
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solve(Scalar timeStep)
{
    grid_->comm().printf("Updating IB positions...\n");
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    solveGammaEqn(timeStep);
    updateProperties(timeStep);
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    computeIbForces(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solveGammaEqn(Scalar timeStep)
{
    auto beta = axi::cicsam::faceInterpolationWeights(u_, gamma_, gradGamma_, timeStep);

    gamma_.savePreviousTimeStep(timeStep, 1.);
    gammaEqn_ = (axi::ddt(gamma_, timeStep) + axi::cicsam::div(u_, gamma_, beta, 0.5, *fluid_) == 0.);

    gammaEqn_.solve();
    grid_->sendMessages(gamma_);
    gamma_.interpolateFaces();
    gradGamma_.computeAxisymmetric(*fluid_);
    grid_->sendMessages(gradGamma_);

    rhoU_.savePreviousTimeStep(timeStep, 1.);
    axi::cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_, beta, rhoU_);
    axi::cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_.oldField(0), beta, rhoU_.oldField(0));
}

void FractionalStepAxisymmetricDFIBMultiphase::updateProperties(Scalar timeStep)
{
    rho_.savePreviousTimeStep(timeStep, 1);

    rho_.computeCells([this](const Cell &c) {
        Scalar g = clamp(gamma_(c), 0., 1.);
        return (1. - g) * rho1_ + g * rho2_;
    });

    rho_.computeFaces([this](const Face &f) {
        Scalar g = clamp(gamma_(f), 0., 1.);
        return (1. - g) * rho1_ + g * rho2_;
    });

    mu_.savePreviousTimeStep(timeStep, 1);

    mu_.computeCells([this](const Cell &c) {
        Scalar g = clamp(gamma_(c), 0., 1.);
        return (1. - g) * mu1_ + g * mu2_;
    });

    mu_.computeFaces([this](const Face &f) {
        Scalar g = clamp(gamma_(f), 0., 1.);
        return (1. - g) * mu1_ + g * mu2_;
    });
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (axi::ddt(rho_, u_, timeStep) + axi::div(rhoU_, u_, 0.) == axi::divSigma(rho_, mu_, p_, u_, 0.));

    uEqn_.solve();
    grid_->sendMessages(u_);

    uEqn_ == axi::divSigma(rho_, mu_, p_, u_, 0.5) - axi::divSigma(rho_, mu_, p_, u_, 0.)
            + ib_->polarVelocityBcs(rho_, u_, u_, timeStep);

    u_.savePreviousIteration();
    uEqn_.solve();

    for(const Cell &c: *fluid_)
    {
        fib_(c) = rho_(c) * (u_(c) - u_.prevIteration()(c)) / timeStep;

        for(const InteriorLink &nb: c.neighbours())
            u_(c) += timeStep / rho_(nb.face()) * p_(nb.face()) * nb.polarOutwardNorm() / c.polarVolume();

        for(const BoundaryLink &bd: c.boundaries())
            u_(c) += timeStep / rho_(bd.face()) * p_(bd.face()) * bd.polarOutwardNorm() / c.polarVolume();

        u_(c) -= timeStep / rho_(c) * p_(c) * Vector2D(c.volume(), 0.) / c.polarVolume();
    }

    grid_->sendMessages(fib_);
    grid_->sendMessages(u_);

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.volumeWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        if(ib_->ibObj(f.lCell().centroid()) || ib_->ibObj(f.rCell().centroid()))
            u_(f) = g * u_(l) + (1. - g) * u_(r);
        else
            u_(f) = g * (u_(l) - timeStep / rho_(l) * (fst_(l) + sg_(l)))
                    + (1. - g) * (u_(r) - timeStep / rho_(r) * (fst_(r) + sg_(r)))
                    + timeStep / rho_(f) * (fst_(f) + sg_(f));
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
                u_(f) = u_(l) - timeStep / rho_(l) * (fst_(l) + sg_(l))
                        + timeStep / rho_(f) * (fst_(f) + sg_(f));
            }
            break;
        case VectorFiniteVolumeField::SYMMETRY:
            for (const Face &f: patch)
            {
                Vector2D tw = f.norm().tangentVec();
                u_(f) = dot(u_(f.lCell()), tw) * tw / tw.magSqr();
            }
            break;
        }
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep / rho_, p_) == axi::src::div(u_));
    pEqn_.solve();
    grid_->sendMessages(p_);
    p_.interpolateFaces();

    gradP_.computeAxisymmetric(rho_, *fluid_);
    grid_->sendMessages(gradP_);
}

void FractionalStepAxisymmetricDFIBMultiphase::correctVelocity(Scalar timeStep)
{
    for(const Face &f: grid_->faces())
        u_(f) -= timeStep / rho_(f) * gradP_(f);

    for(const Cell &c: *fluid_)
        u_(c) -= timeStep / rho_(c) * gradP_(c);

    grid_->sendMessages(u_);
}

void FractionalStepAxisymmetricDFIBMultiphase::computeIbForces(Scalar timeStep)
{
    for(auto &ibObj: *ib_)
    {
        Vector2D fh(0., 0.);

        for(const Cell &c: ibObj->solidCells())
            fh -= fib_(c) * c.polarVolume() * 2. * M_PI;

        for(const Cell &c: ibObj->ibCells())
            fh -= fib_(c) * c.polarVolume() * 2. * M_PI;

        fh = grid_->comm().sum(fh);

        Vector2D fw(0., 0.);

        if(ibObj->shape().type() == Shape2D::CIRCLE)
        {
            //- Assume spherical
            const Circle &circ = static_cast<const Circle&>(ibObj->shape());

            if(circ.centroid().x == 0.)
            {
                Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);
                fw = (ibObj->rho - FractionalStep::rho_) * vol * g_;
            }
        }

        ibObj->applyForce(fh + fw);
    }
}
