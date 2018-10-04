#include "Math/Algorithm.h"
#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricExplicitDivergence.h"
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
      gammaSrc_(*addField<Scalar>("gammaSrc", fluid_)),
      rho_(*addField<Scalar>("rho", fluid_)),
      mu_(*addField<Scalar>("mu", fluid_)),
      rhoU_(*addField<Vector2D>("rhoU", fluid_)),
      sg_(*addField<Vector2D>("sg", fluid_)),
      gradGamma_(static_cast<ScalarGradient&>(*addField<Vector2D>(std::make_shared<ScalarGradient>(gamma_, fluid_)))),
      gradRho_(static_cast<ScalarGradient&>(*addField<Vector2D>(std::make_shared<ScalarGradient>(rho_, fluid_)))),
      fst_(input, grid, fluid_, ib_),
      gammaEqn_(input, gamma_, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", FractionalStep::rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", FractionalStep::rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", FractionalStep::mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", FractionalStep::mu_);

    //- Set axisymmetric
    fst_.setAxisymmetric(true);

    //- Register surface tension fields
    addField<Scalar>(fst_.gammaTilde());
    addField<Scalar>(fst_.kappa());
    addField<Vector2D>(fst_.gradGammaTilde());
    addField<Vector2D>(fst_.n());
    addField<Vector2D>(fst_.fst());
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
    gammaEqn_ = (axi::ddt(gamma_, timeStep) + axi::cicsam::div(u_, gamma_, beta, 0., *fluid_) == 0.);
    gammaEqn_.solve();
    gamma_.sendMessages();

    gamma_.savePreviousIteration();
    fst_.computeContactLineExtension(gamma_);

    for(const Cell &c: *fluid_)
        gammaSrc_(c) = (gamma_(c) - gamma_.prevIteration()(c)) / timeStep;

    gammaEqn_ == axi::src::src(gammaSrc_);

    gammaEqn_.solve();
    gamma_.sendMessages();
    gamma_.interpolateFaces();

    rhoU_.savePreviousTimeStep(timeStep, 2.);
    axi::cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_, beta, rhoU_.oldField(0));
    axi::cicsam::computeMomentumFlux(rho1_, rho2_, u_, gamma_.oldField(0), beta, rhoU_.oldField(1));

    gradGamma_.computeAxisymmetric(*fluid_);
    gradGamma_.sendMessages();
}

void FractionalStepAxisymmetricDFIBMultiphase::updateProperties(Scalar timeStep)
{
    //- Update rho and mu
    rho_.savePreviousTimeStep(timeStep, 1);

    rho_.computeCells([this](const Cell &c) {
        return rho1_ + clamp(gamma_(c), 0., 1.) * (rho2_ - rho1_);
    });

    rho_.computeFaces([this](const Face &f) {
        return rho1_ + clamp(gamma_(f), 0., 1.) * (rho2_ - rho1_);
    });

    mu_.savePreviousTimeStep(timeStep, 1);

    mu_.computeCells([this](const Cell &c) {
        return rho_(c) / (rho1_ / mu1_ + clamp(gamma_(c), 0., 1.) * (rho2_ / mu2_ - rho1_ / mu1_));
    });

    mu_.computeFaces([this](const Face &f) {
        return rho_(f) / (rho1_ / mu1_ + clamp(gamma_(f), 0., 1.) * (rho2_ / mu2_ - rho1_ / mu1_));
    });

    //- Update the gravitational source term
    //- Update the gravitational source term
    gradRho_.computeFaces();

    sg_.computeFaces([this](const Face &face) {
        return dot(g_, -face.centroid()) * gradRho_(face);
    });

    sg_.faceToCellAxisymmetric(rho_, rho_, *fluid_);
    sg_.sendMessages();

    //- Update the surface tension force
    fst_.computeFaceInterfaceForces(gamma_, gradGamma_);
    fst_.fst()->faceToCellAxisymmetric(rho_, rho_, *fluid_);
    fst_.fst()->fill(Vector2D(0., 0.), ib_->localSolidCells());
    fst_.fst()->sendMessages();
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solveUEqn(Scalar timeStep)
{
    gradP_.computeAxisymmetric(rho_, rho_.oldField(0), *fluid_);

    const VectorFiniteVolumeField &fst = *fst_.fst();

    u_.savePreviousTimeStep(timeStep, 2);
    uEqn_ = (axi::ddt(rho_, u_, timeStep) + axi::dive(rhoU_, u_, 0.5)
             == axi::laplacian(mu_, u_, 0.) + axi::src::src(sg_ + fst - gradP_));

    uEqn_.solve();
    u_.sendMessages();

    uEqn_ == axi::laplacian(mu_, u_, 0.5) - axi::laplacian(mu_, u_, 0.)
            + ib_->polarVelocityBcs(rho_, u_, u_, timeStep);

    u_.savePreviousIteration();
    uEqn_.solve();

    for(const Cell &c: *fluid_)
    {
        fib_(c) = rho_(c) * (u_(c) - u_.prevIteration()(c)) / timeStep;
        u_(c) += timeStep / rho_(c) * gradP_(c);
    }

    fib_.sendMessages();
    u_.sendMessages();

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.distanceWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        if(ib_->ibObj(f.lCell().centroid()) || ib_->ibObj(f.rCell().centroid()))
            u_(f) = g * (u_(l) - timeStep / rho_(l) * sg_(l))
                    + (1. - g) * (u_(r) - timeStep / rho_(r) * sg_(r))
                    + timeStep / rho_(f) * sg_(f);
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
    p_.sendMessages();
    p_.setBoundaryFaces();

    gradP_.computeAxisymmetric(rho_, rho_, *fluid_);
    gradP_.sendMessages();
}

void FractionalStepAxisymmetricDFIBMultiphase::correctVelocity(Scalar timeStep)
{
    for(const Face &f: grid_->faces())
        u_(f) -= timeStep / rho_(f) * gradP_(f);

    for(const Cell &c: *fluid_)
        u_(c) -= timeStep / rho_(c) * gradP_(c);

    u_.sendMessages();
}

void FractionalStepAxisymmetricDFIBMultiphase::computeIbForces(Scalar timeStep)
{
    struct ContactLine
    {
        ContactLine() = default;

        ContactLine(const Point2D &pt, Scalar gamma, const Vector2D &n = Vector2D(0., 0.), const Vector2D &t = Vector2D(0., 0.))
            :
              pt(pt), gamma(gamma), n(n), t(t)
        {}

        Point2D pt;

        Scalar gamma;

        Vector2D n, t;
    };

    std::vector<ContactLine> tmp;

    for(auto &ibObj: *ib_)
    {
        Vector2D fh(0., 0.);

        for(const Cell &c: ibObj->cells())
        {
//            Vector2D DuDt = (rho_(c) * u_(c) - rho_.oldField(0)(c) * u_.oldField(0)(c)) / timeStep;

//            for(const InteriorLink &nb: c.neighbours())
//            {
//                const Face &f = nb.face();
//                DuDt += 0.5 * dot(outer(rhoU_.oldField(0)(f), u_.oldField(0)(f)), nb.polarOutwardNorm()) / c.polarVolume();
//                DuDt += 0.5 * dot(outer(rhoU_.oldField(1)(f), u_.oldField(1)(f)), nb.polarOutwardNorm()) / c.polarVolume();
//            }

//            for(const BoundaryLink &bd: c.boundaries())
//            {
//                const Face &f = bd.face();
//                DuDt += 0.5 * dot(outer(rhoU_.oldField(0)(f), u_.oldField(0)(f)), bd.polarOutwardNorm()) / c.polarVolume();
//                DuDt += 0.5 * dot(outer(rhoU_.oldField(1)(f), u_.oldField(1)(f)), bd.polarOutwardNorm()) / c.polarVolume();
//            }

            fh -= fib_(c) * 2. * M_PI * c.polarVolume();
        }

        fh = grid_->comm().sum(fh);

        Vector2D fc(0., 0.), fb(0., 0.), fw(0., 0.);

        if(ibObj->shape().type() == Shape2D::CIRCLE && ibObj->shape().centroid().x == 0.)
        {
            //- Assume spherical
            const Circle &circ = static_cast<const Circle&>(ibObj->shape());

            tmp.clear();
            tmp.reserve(ibObj->solidCells().size() + 2);

            for(const Cell &c: ibObj->solidCells())
            {
                bool computeCl = false;

                for(const CellLink &nb: c.neighbours())
                    if(!ib_->ibObj(nb.cell().centroid()))
                        computeCl = true;

                if(!computeCl)
                    continue;

                auto st = fst_.contactLineStencil(c.centroid(), gamma_);

                tmp.emplace_back(st.cl()[1], st.gamma(), st.ncl(), st.tcl());
            }

            tmp = grid_->comm().allGatherv(tmp);

            std::sort(tmp.begin(), tmp.end(), [&ibObj](const ContactLine &lhs, const ContactLine &rhs)
            { return (lhs.pt - ibObj->shape().centroid()).angle() < (rhs.pt - ibObj->shape().centroid()).angle(); });

            auto buoyancyForce = [this, &circ](const ContactLine &cl1, const ContactLine &cl2)
            {
                Scalar g1 = cl1.gamma;
                Scalar g2 = cl2.gamma;

                Scalar th1 = (cl1.pt - circ.centroid()).angle();
                Scalar th2 = (cl2.pt - circ.centroid()).angle();
                th2 = th2 < th1 ? th2 + 2. * M_PI : th2;

                Scalar rgh1 = (rho1_ + clamp(g1, 0., 1.) * (rho2_ - rho1_)) * dot(cl1.pt, g_);
                Scalar rgh2 = (rho1_ + clamp(g2, 0., 1.) * (rho2_ - rho1_)) * dot(cl2.pt, g_);

                //- Find the buoyancy
                Scalar r = circ.radius();
                return -Vector2D(0.,
                                 M_PI * r * r * (-rgh1 * th1 * std::cos(2. * th1) + rgh1 * th2 * std::cos(2. * th1)
                                                 + rgh1 * std::sin(2. * th1) / 2. - rgh1 * std::sin(2. * th2) / 2.
                                                 + rgh2 * th1 * std::cos(2. * th2) - rgh2 * th2 * std::cos(2. * th2)
                                                 - rgh2 * std::sin(2. * th1) / 2. + rgh2 * std::sin(2. * th2) / 2.) / (2. * (th1 - th2))
                                 );
            };

            for(int i = 0; i < tmp.size() - 1; ++i)
            {
                const auto &cl1 = tmp[i];
                const auto &cl2 = tmp[(i + 1) % tmp.size()];

                if(i == 0)
                    fb += buoyancyForce(ContactLine(ibObj->nearestIntersect(Point2D(0., cl1.pt.y)), cl1.gamma), cl1);

                if(i == tmp.size() - 2)
                    fb += buoyancyForce(cl2, ContactLine(ibObj->nearestIntersect(Point2D(0., cl2.pt.y)), cl2.gamma));

                fb += buoyancyForce(cl1, cl2);

                Scalar g1 = cl1.gamma;
                Scalar g2 = cl2.gamma;
                Scalar th1 = (cl1.pt - circ.centroid()).angle();
                Scalar th2 = (cl2.pt - circ.centroid()).angle();

                // check if contact line exists between two points
                if(g1 < 0.5 == g2 <= 0.5)
                    continue;

                //- Interpolate along an arc
                Scalar alpha = (0.5 - g1) / (g2 - g1);
                Scalar th = th1 + alpha * (th2 - th1);
                Point2D pt = circ.centroid() + (cl1.pt - circ.centroid()).rotate(th - th1);

                // construct new contact line at point to find capillary force
                auto cl = CelesteImmersedBoundary::ContactLineStencil(*ibObj, pt, fst_.theta(*ibObj), gamma_);

                if(cl.isValid())
                    fc += 2. * M_PI * pt.x * fst_.sigma() * cl.tcl();
            }

            if(circ.centroid().x == 0.)
            {
                Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);
                fw = ibObj->rho * vol * g_;
            }

            if(grid_->comm().isMainProc())
            {
                Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);

                std::cout << "Hydro force = " << fh << "\n"
                          << "Weight = " << fw + FractionalStep::rho_ * vol * g_ << "\n"
                          << "Capillary force = " << fc << "\n"
                          << "Buoyancy force = " << fb << " == " << FractionalStep::rho_ * vol * g_ << "\n";
            }
        }

        ibObj->applyForce(fh + fw + fc + fb);
    }
}
