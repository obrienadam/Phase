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

    grid_->comm().printf("Solving gamma equation and updating properties...\n");
    solveGammaEqn(timeStep);
    updateProperties(timeStep);

    grid_->comm().printf("Solving momentum and pressure equation and correting velocity...\n");
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Computing IB forces...\n");
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

    Scalar error = gammaEqn_.solve();
    gamma_.sendMessages();
    gamma_.interpolateFaces();

    gradGamma_.computeAxisymmetric(*fluid_);
    gradGamma_.sendMessages();

    return error;
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
    gradRho_.computeFaces();

    sg_.computeFaces([this](const Face &face) {
        return dot(g_, -face.centroid()) * gradRho_(face);
    });

    sg_.faceToCellAxisymmetric(rho_, rho_, *fluid_);
    sg_.sendMessages();

    //- Update the surface tension
    fst_.computeFaceInterfaceForces(gamma_, gradGamma_);
    fst_.fst()->faceToCellAxisymmetric(rho_, rho_, *fluid_);
    fst_.fst()->sendMessages();
}

Scalar FractionalStepAxisymmetricDFIBMultiphase::solveUEqn(Scalar timeStep)
{
    gradP_.computeAxisymmetric(rho_, rho_.oldField(0), *fluid_);
    //gradP_.fill(Vector2D(0., 0.), ib_->localSolidCells());
    //gradP_.fill(Vector2D(0., 0.), ib_->localIbCells());
    gradP_.sendMessages();

    const VectorFiniteVolumeField &fst = *fst_.fst();

    u_.savePreviousTimeStep(timeStep, 2);
    uEqn_ = (rho_ * axi::ddt(u_, timeStep) + rho_ * axi::dive(u_, u_, 0.5)
             == axi::laplacian(mu_, u_, 0.5) + axi::src::src(sg_ + fst - gradP_));

    Scalar error = uEqn_.solve();
    u_.sendMessages();

    fibEqn_ = ib_->computeForcingTerm(u_, timeStep, fib_);
    fibEqn_.solve();
    fib_.sendMessages();

    for(const Cell &c: *fluid_)
        u_(c) += timeStep * (fib_(c) + gradP_(c) / rho_(c));

    u_.sendMessages();

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.distanceWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        if(ib_->ibObj(l) || ib_->ibObj(r))
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

Scalar FractionalStepAxisymmetricDFIBMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep / rho_, p_) == axi::src::div(u_));

    Scalar error = pEqn_.solve();
    p_.sendMessages();
    p_.setBoundaryFaces();

    gradP_.computeAxisymmetric(rho_, rho_, *fluid_);
    gradP_.sendMessages();

    return error;
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
    for(auto &ibObj: *ib_)
    {
        if(ibObj->shape().type() != Shape2D::CIRCLE || ibObj->shape().centroid().x != 0.)
            throw Exception("FractionalStepAxisymmetricDFIBMultiphase", "computeIbForces", "oncly circles centered at r = 0 supported.");

        Vector2D fh(0., 0.);
        for(const Cell &c: ibObj->cells())
        {
            fh += rho_(c) * (u_(c) - u_.oldField(0)(c)) * c.polarVolume() / timeStep;

            for(const InteriorLink &nb: c.neighbours())
            {
                Scalar flux0 = rho_(c) * dot(u_.oldField(0)(nb.face()), nb.polarOutwardNorm()) / 2.;
                Scalar flux1 = rho_(c) * dot(u_.oldField(1)(nb.face()), nb.polarOutwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(nb.cell())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(nb.cell());
            }

            for(const BoundaryLink &bd: c.boundaries())
            {
                Scalar flux0 = rho_(c) * dot(u_.oldField(0)(bd.face()), bd.polarOutwardNorm()) / 2.;
                Scalar flux1 = rho_(c) * dot(u_.oldField(1)(bd.face()), bd.polarOutwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(bd.face())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(bd.face());
            }

            fh -= rho_(c) * (fib_(c) + (*fst_.fst())(c) / rho_(c) + g_) * c.polarVolume();
        }

        fh = grid_->comm().sum(2. * M_PI * fh);

        Vector2D fc(0., 0.), fw(0., 0.);

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

            auto st = CelesteAxisymmetricImmersedBoundary::ContactLineStencil(*ibObj, c.centroid(), fst_.theta(*ibObj), gamma_);
            Scalar beta = (st.cl()[1] - ibObj->shape().centroid()).angle();
            contactLines_.emplace_back(ContactLine{st.cl()[1], beta, st.gamma(), st.ncl(), st.tcl()});
        }

        contactLines_ = grid_->comm().allGatherv(contactLines_);

        std::sort(contactLines_.begin(), contactLines_.end(), [](const ContactLine &lhs, const ContactLine &rhs)
        { return lhs.beta < rhs.beta; });

        //- Assume spherical
        const Circle &circ = static_cast<const Circle&>(ibObj->shape());

        //- Assumes the sphere is centered on the axis
        for(size_t i = 0; i < contactLines_.size() - 1; ++i)
        {
            const auto &cl1 = contactLines_[i];
            const auto &cl2 = contactLines_[(i + 1) % contactLines_.size()];

            Scalar g1 = cl1.gamma;
            Scalar g2 = cl2.gamma;
            Scalar beta1 = cl1.beta;
            Scalar beta2 = cl2.beta < cl1.beta ? cl2.beta + 2. * M_PI : cl2.beta;


            // check if contact line exists between two points
            //- sharp method
            if((g1 < 0.5) != (g2 <= 0.5))
            {
                Scalar alpha = (0.5 - g2) / (g1 - g2);
                Scalar beta = alpha * beta1 + (1. - alpha) * beta2;
                Vector2D tcl = alpha < 0.5 ? cl1.tcl.rotate(beta - beta1) : cl2.tcl.rotate(beta - beta2);

                if(grid_->comm().isMainProc())
                    std::cout << "Beta = " << beta * 180. / M_PI << "\n"
                              << "Contact line = " << tcl << "\n";

                Point2D pt = circ.centroid() + (cl1.pt - circ.centroid()).rotate(beta - beta1);
                Scalar r = pt.x;

                fc += fst_.sigma() * 2. * M_PI * r * tcl;
            }
        }

        Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);

        fw = ibObj->rho * vol * g_;
        if(grid_->comm().isMainProc())
        {
            std::cout << "Hydro force = " << fh << "\n"
                      << "Weight = " << fw << "\n"
                      << "Capillary force = " << fc << "\n"
                      << "Net force = " << fh + fw + fc << "\n";
        }

        ibObj->applyForce((fh + fw + fc) * ibObj->mass() / (ibObj->rho * vol));
    }
}

void FractionalStepAxisymmetricDFIBMultiphase::computeFieldExtenstions(Scalar timeStep)
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

                CelesteImmersedBoundary::ContactLineStencil cl(*ibObj, c.centroid(), fst_.theta(*ibObj), gamma_);
                CelesteImmersedBoundary::ContactLineStencil stn(*ibObj, c.centroid(), M_PI_2, gamma_);


                Scalar ubn = dot(ibObj->velocity(bp), ns);
                Scalar abn = dot(ibObj->acceleration(bp), ns);
                Scalar rhob = cl.interpolate(rho_);
                Scalar dRho = (rho_(c) - rhob) / (c.centroid() - bp).mag();
                Scalar pb = stn.interpolate(p_);
                Scalar dP = -(2 * ubn * ubn * dRho + rhob * abn);

                if(!std::isnan(dP))
                {
                    p_(c) = pb + dP * (c.centroid() - bp).mag();
                }
            }
        }

    gradP_.computeAxisymmetric(rho_, rho_, *fluid_);
    gradP_.sendMessages();
}
