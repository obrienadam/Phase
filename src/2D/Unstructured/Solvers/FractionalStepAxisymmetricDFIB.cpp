#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"
#include "FiniteVolume/Motion/SolidBodyMotion.h"
#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricStressTensor.h"
#include "FiniteVolume/Discretization/AxisymmetricSource.h"

#include "FractionalStepAxisymmetricDFIB.h"

FractionalStepAxisymmetricDFIB::FractionalStepAxisymmetricDFIB(const Input &input,
                                                               const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepAxisymmetric(input, grid),
      fib_(*addField<Vector2D>("fb", fluid_)),
      fibEqn_(input, fib_, "fbEqn")
{
    ib_ = std::make_shared<DirectForcingImmersedBoundary>(input, grid, fluid_);
    addField<int>(ib_->cellStatus());

    for(auto &ibObj: *ib_)
    {
        auto motion = std::dynamic_pointer_cast<SolidBodyMotion>(ibObj->motion());

        if(motion)
            motion->setMotionConstraint(Vector2D(0., 1.));
    }
}

Scalar FractionalStepAxisymmetricDFIB::solve(Scalar timeStep)
{
    grid_->comm().printf("Updating IB positions...\n");
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    computeIbForces(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0.;
}

std::shared_ptr<const ImmersedBoundary> FractionalStepAxisymmetricDFIB::ib() const
{
    return ib_;
}

Scalar FractionalStepAxisymmetricDFIB::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (axi::ddt(u_, timeStep) + axi::div(u_, u_, 0.)
             == axi::divSigma(mu_ / rho_, p_, u_, 0.));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u_);

    uEqn_ == axi::divSigma(mu_ / rho_, p_, u_, 0.5) - axi::divSigma(mu_ / rho_, p_, u_, 0.)
            + ib_->polarVelocityBcs(u_, u_, timeStep);

    u_.savePreviousIteration();
    uEqn_.solve();

    for(const Cell &c: *fluid_)
        fib_(c) = (u_(c) - u_.prevIteration()(c)) / timeStep;

    grid_->sendMessages(fib_);

    for(const Cell &c: *fluid_)
        u_(c) += timeStep * gradP_(c);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    return error;
}

void FractionalStepAxisymmetricDFIB::computeIbForces(Scalar timeStep)
{
    for(auto &ibObj: *ib_)
    {
        Vector2D fh(0., 0.);

        for(const Cell &c: ibObj->solidCells())
            fh -= rho_ * fib_(c) * c.polarVolume() * 2. * M_PI;

        for(const Cell &c: ibObj->ibCells())
            fh -= rho_ * fib_(c) * c.polarVolume() * 2. * M_PI;

        fh = grid_->comm().sum(fh);

        Vector2D fw(0., 0.);

        if(ibObj->shape().type() == Shape2D::CIRCLE)
        {
            //- Assume spherical
            const Circle &circ = static_cast<const Circle&>(ibObj->shape());

            if(circ.centroid().x == 0.)
            {
                Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);
                fw = (ibObj->rho - rho_) * vol * g_;
            }


            std::cout << "Cd = " << 2. * fh.y / (rho_ * M_PI * std::pow(circ.radius(), 2)) << "\n";
        }

        ibObj->applyForce(fh + fw);
    }
}
