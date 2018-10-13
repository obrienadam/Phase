#include "Math/Algorithm.h"

#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"
#include "FiniteVolume/Multiphase/CelesteImmersedBoundary.h"
#include "FiniteVolume/Discretization/TimeDerivative.h"
#include "FiniteVolume/Discretization/Divergence.h"
#include "FiniteVolume/Discretization/ExplicitDivergence.h"
#include "FiniteVolume/Discretization/Laplacian.h"
#include "FiniteVolume/Discretization/Source.h"
#include "FiniteVolume/Discretization/Cicsam.h"

#include "FractionalStepDFIBMultiphase.h"

FractionalStepDirectForcingMultiphase::FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepDFIB(input, grid),
      gamma_(*addField<Scalar>(input, "gamma", fluid_)),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      gammaSrc_(*addField<Scalar>("gammaSrc", fluid_)),
      sg_(*addField<Vector2D>("sg", fluid_)),
      gradGamma_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(gamma_, fluid_)))),
      gradRho_(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(rho_, fluid_)))),
      fst_(std::make_shared<CelesteImmersedBoundary>(input, grid_, fluid_, ib_)),
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
        capillaryTimeStep_ = std::min(
                    capillaryTimeStep_,
                    std::sqrt(
                        (rho1_ + rho2_) * std::pow(delta, 3) / (4. * M_PI * fst_->sigma())
                        ));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);

    addField(fst_->fst());
    addField(fst_->kappa());
    addField(fst_->gammaTilde());
    addField<Vector2D>(fst_->gradGammaTilde());
    addField(fst_->n());
}

void FractionalStepDirectForcingMultiphase::initialize()
{
    FractionalStepDFIB::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma_.compute(*fluid_);
    updateProperties(0.);
}

Scalar FractionalStepDirectForcingMultiphase::solve(Scalar timeStep)
{
    //- Perform field extension
    grid_->comm().printf("Updating IB positions and cell categories...\n");
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    grid_->comm().printf("Solving gamma equation...\n");
    solveGammaEqn(timeStep);

    grid_->comm().printf("Updating physical properties...\n");
    updateProperties(timeStep);

    grid_->comm().printf("Solving momentum equation...\n");
    solveUEqn(timeStep);

    grid_->comm().printf("Solving pressure equation and correcting velocities...\n");
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Computing IB forces...\n");
    computeIbForces(timeStep);
    ib_->applyCollisionForce(true);

    //    grid_->comm().printf("Performing field extensions...\n");
    //    computeFieldExtenstions(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepDirectForcingMultiphase::solveGammaEqn(Scalar timeStep)
{
    auto beta = cicsam::faceInterpolationWeights(u_, gamma_, gradGamma_, timeStep);

    //- Predictor
    gamma_.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma_, timeStep) + cicsam::div(u_, gamma_, beta, 0.) == 0.);

    Scalar error = gammaEqn_.solve();
    gamma_.sendMessages();

    //- Corrector
    gamma_.savePreviousIteration();
    fst_->computeContactLineExtension(gamma_);

    for(const Cell &c: *fluid_)
        gammaSrc_(c) = (gamma_(c) - gamma_.prevIteration()(c)) / timeStep;

    gammaEqn_ == cicsam::div(u_, gamma_, beta, 0.5) - cicsam::div(u_, gamma_, beta, 0.) + src::src(gammaSrc_);

    error = gammaEqn_.solve();
    gamma_.sendMessages();
    gamma_.interpolateFaces();

    //- Update the gradient
    gradGamma_.compute(*fluid_);
    gradGamma_.sendMessages();

    return error;
}

Scalar FractionalStepDirectForcingMultiphase::solveUEqn(Scalar timeStep)
{
    const auto &fst = *fst_->fst();
    gradP_.faceToCell(rho_, rho_.oldField(0), *fluid_);

    //- Explicit predictor
    u_.savePreviousTimeStep(timeStep, 2);
    uEqn_ = (rho_ * fv::ddt(u_, timeStep) + rho_ * fv::dive(u_, u_, 0.5)
             == fv::laplacian(mu_, u_, 0.) + src::src(fst + sg_ - gradP_));

    Scalar error = uEqn_.solve();
    u_.sendMessages();

    //- Semi-implicit corrector
    u_.savePreviousIteration();
    uEqn_ == fv::laplacian(mu_, u_, 0.5) - fv::laplacian(mu_, u_, 0.)
            + ib_->velocityBcs(rho_, u_, u_, timeStep);
    error = uEqn_.solve();

    for(const Cell& c: *fluid_)
    {
        fb_(c) = rho_(c) * (u_(c) - u_.prevIteration()(c)) / timeStep;
        u_(c) += timeStep / rho_(c) * gradP_(c);
    }

    fb_.sendMessages();
    u_.sendMessages();

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.volumeWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        if(ib_->ibObj(f.lCell().centroid()) || ib_->ibObj(f.rCell().centroid()))
            u_(f) = g * u_(l) + (1. - g) * u_(r);
        else
            u_(f) = g * (u_(l) - timeStep / rho_(l) * (fst(l) + sg_(l)))
                    + (1. - g) * (u_(r) - timeStep / rho_(r) * (fst(r) + sg_(r)))
                    + timeStep / rho_(f) * (fst(f) + sg_(f));
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
                u_(f) = u_(f.lCell()).tangentialComponent(f.norm());
            break;
        }

    return error;
}

Scalar FractionalStepDirectForcingMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p_) == src::div(u_));

    Scalar error = pEqn_.solve();
    p_.sendMessages();
    p_.setBoundaryFaces();

    gradP_.computeFaces();
    gradP_.faceToCell(rho_, rho_, *fluid_);
    gradP_.sendMessages();

    return error;
}

void FractionalStepDirectForcingMultiphase::updateProperties(Scalar timeStep)
{
    //- Update density
    rho_.savePreviousTimeStep(timeStep, 1);

    rho_.computeCells([this](const Cell &c) {
        return rho1_ + clamp(gamma_(c), 0., 1.) * (rho2_ - rho1_);
    });

    rho_.computeFaces([this](const Face &f) {
        return rho1_ + clamp(gamma_(f), 0., 1.) * (rho2_ - rho1_);
    });

    //- Update the gravitational source term
    gradRho_.computeFaces();

    for (const Face &face: grid_->faces())
        sg_(face) = -dot(g_, face.centroid()) * gradRho_(face);

    sg_.faceToCell(rho_, rho_, *fluid_);
    sg_.fill(Vector2D(0., 0.), ib_->localSolidCells());
    sg_.sendMessages();

    //- Update viscosity from kinematic viscosity
    mu_.savePreviousTimeStep(timeStep, 1);

    mu_.computeCells([this](const Cell &c) {
        return rho_(c) / (rho1_ / mu1_ + clamp(gamma_(c), 0., 1.) * (rho2_ / mu2_ - rho1_ / mu1_));
    });

    mu_.computeFaces([this](const Face &f) {
        return rho_(f) / (rho1_ / mu1_ + clamp(gamma_(f), 0., 1.) * (rho2_ / mu2_ - rho1_ / mu1_));
    });

    //- Update the surface tension
    fst_->computeFaceInterfaceForces(gamma_, gradGamma_);
    fst_->fst()->faceToCell(rho_, rho_, *fluid_);
    fst_->fst()->fill(Vector2D(0., 0.), ib_->localSolidCells());
    fst_->fst()->sendMessages();
}

void FractionalStepDirectForcingMultiphase::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: *fluid_)
        u_(cell) -= timeStep / rho_(cell) * gradP_(cell);

    u_.sendMessages();

    for (const Face &face: grid_->faces())
        u_(face) -= timeStep / rho_(face) * gradP_(face);
}

void FractionalStepDirectForcingMultiphase::computeIbForces(Scalar timeStep)
{
    for(auto &ibObj: *ib_)
    {
        contactLines_.clear();

        auto computeStress = [&ibObj](const Cell &c)
        {
            for(const CellLink &nb: c.neighbours())
                if(!ibObj->isInIb(nb.cell().centroid()))
                    return true;

            for(const CellLink &nb: c.diagonals())
                if(!ibObj->isInIb(nb.cell().centroid()))
                    return true;

            return false;
        };

        for(const Cell &c: ibObj->solidCells())
        {
            if(!computeStress(c))
                continue;

            Point2D pt = ibObj->nearestIntersect(c.centroid());
            Scalar beta = (pt - ibObj->shape().centroid()).angle();

            CelesteImmersedBoundary::ContactLineStencil st1(*ibObj, pt, fst_->theta(*ibObj), gamma_);

            Scalar rho = st1.interpolate(rho_);
            Scalar rgh = rho * dot(g_, pt - ibObj->shape().centroid());
            Scalar gamma = st1.gamma();
            Vector2D ncl = st1.ncl();
            Vector2D tcl = st1.tcl();

            contactLines_.push_back(ContactLine{pt, beta, rho, rgh, gamma, ncl, tcl});
        }

        contactLines_ = grid_->comm().allGatherv(contactLines_);

        std::sort(contactLines_.begin(), contactLines_.end(), [](const ContactLine &lhs, const ContactLine &rhs)
        { return lhs.beta < rhs.beta; });

        Vector2D fb(0., 0.), fc(0., 0.);

        if(ibObj->shape().type() == Shape2D::CIRCLE)
        {
            const Circle &circ = static_cast<const Circle&>(ibObj->shape());
            Scalar r = circ.radius();

            for(auto i = 0; i < contactLines_.size(); ++i)
            {
                const ContactLine &stA = contactLines_[i];
                const ContactLine &stB = contactLines_[(i + 1) % contactLines_.size()];

                Scalar tA = stA.beta;
                Scalar tB = stB.beta < tA ? stB.beta + 2. * M_PI : stB.beta;

                Scalar pA = stA.rgh;
                Scalar pB = stB.rgh;

                fb += Vector2D(
                            r*(pA*tA*sin(tA) - pA*tB*sin(tA) + pA*cos(tA) - pA*cos(tB) - pB*tA*sin(tB) + pB*tB*sin(tB) - pB*cos(tA) + pB*cos(tB))/(tA - tB),
                            r*(-pA*tA*cos(tA) + pA*tB*cos(tA) + pA*sin(tA) - pA*sin(tB) + pB*tA*cos(tB) - pB*tB*cos(tB) - pB*sin(tA) + pB*sin(tB))/(tA - tB)
                            );

                Scalar gA = stA.gamma;
                Scalar gB = stB.gamma;
                Scalar sigma = fst_->sigma();

                //- sharp method
                if(gA <= 0.5 != gB < 0.5)
                {
                    Scalar alpha = (0.5 - gB) / (gA - gB);
                    Scalar t = alpha * tA + (1. - alpha) * tB;
                    Vector2D tcl = alpha < 0.5 ? stA.tcl.rotate(t - tA) : stB.tcl.rotate(t - tB);
                    fc += sigma * tcl;
                }
            }
        }

        //- Compute the hydro force from the ib force
        Vector2D fh(0., 0.);
        for(const Cell &c: ibObj->cells())
        {
            fh += rho_(c) * (u_(c) - u_.oldField(0)(c)) * c.volume() / timeStep;

            for(const InteriorLink &nb: c.neighbours())
            {
                Scalar flux0 = rho_(c) * dot(u_.oldField(0)(nb.face()), nb.outwardNorm()) / 2.;
                Scalar flux1 = rho_(c) * dot(u_.oldField(1)(nb.face()), nb.outwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(nb.cell())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(nb.cell());
            }

            for(const BoundaryLink &bd: c.boundaries())
            {
                Scalar flux0 = rho_(c) * dot(u_.oldField(0)(bd.face()), bd.outwardNorm()) / 2.;
                Scalar flux1 = rho_(c) * dot(u_.oldField(1)(bd.face()), bd.outwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(bd.face())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(bd.face());
            }

            fh -= fb_(c) * c.volume();
        }

        fh = grid_->comm().sum(fh);

        Vector2D fw = ibObj->rho * ibObj->shape().area() * g_;

        if(grid_->comm().isMainProc())
            std::cout << "Buoyancy force = " << fb << "\n"
                      << "Hydrodynamic force = " << fh << "\n"
                      << "Capillary force = " << fc << "\n"
                      << "Weight = " << fw << "\n"
                      << "Net = " << fh + fb + fc + fw << "\n";


        ibObj->applyForce(fh + fb + fc + fw);
    }
}

void FractionalStepDirectForcingMultiphase::computeFieldExtenstions(Scalar timeStep)
{
    auto extend = [](const ImmersedBoundaryObject &ibObj, const Cell &c)
    {
        for(const CellLink &nb: c.neighbours())
            if(!ibObj.isInIb(nb.cell().centroid()))
                return true;

        for(const CellLink &nb: c.diagonals())
            if(!ibObj.isInIb(nb.cell().centroid()))
                return true;

        return false;
    };

    for(const auto &ibObj: *ib_)
        for(const Cell &c: ibObj->solidCells())
        {
            if(extend(*ibObj, c))
            {
                Point2D bp = ibObj->nearestIntersect(c.centroid());
                Vector2D ns = ibObj->nearestEdgeUnitNormal(c.centroid());

                CelesteImmersedBoundary::ContactLineStencil cl(*ibObj, c.centroid(), fst_->theta(*ibObj), gamma_);
                CelesteImmersedBoundary::ContactLineStencil stn(*ibObj, c.centroid(), M_PI_2, gamma_);


                Scalar ubn = dot(ibObj->velocity(bp), ns);
                Scalar abn = dot(ibObj->acceleration(bp), ns);
                Scalar rhob = cl.interpolate(rho_);
                Scalar dRho = (rho_(c) - rhob) / (c.centroid() - bp).mag();

                // Convert to static pressure first
                Scalar pb = stn.interpolate(p_);
                Scalar dP = -(2 * ubn * ubn * dRho + rhob * abn);

                if(!std::isnan(dP))
                {
                    p_(c) = pb + dP * (c.centroid() - bp).mag();
                }
            }
        }

    gradP_.computeFaces();
    gradP_.faceToCell(rho_, rho_, *fluid_);
    gradP_.sendMessages();
}
