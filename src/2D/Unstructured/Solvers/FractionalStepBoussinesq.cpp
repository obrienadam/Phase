#include "FiniteVolume/Discretization/TimeDerivative.h"
#include "FiniteVolume/Discretization/Divergence.h"
#include "FiniteVolume/Discretization/Laplacian.h"
#include "FiniteVolume/Discretization/Source.h"

#include "FractionalStepBoussinesq.h"

FractionalStepBoussinesq::FractionalStepBoussinesq(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStep(input, grid),
      T(*addField<Scalar>(input, "T")),
      TEqn_(input, T, "TEqn")
{
    alpha_ = input.caseInput().get<Scalar>("Properties.alpha", 0.00369);
    T0_ = input.caseInput().get<Scalar>("Properties.T0", 273.);
    kappa_ = input.caseInput().get<Scalar>("Properties.kappa", 0.02);
}

Scalar FractionalStepBoussinesq::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveTEqn(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepBoussinesq::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u_, timeStep) + fv::div(u_, u_, 0.5)
             == fv::laplacian(mu_ / rho_, u_, 0.5) - src::src(gradP_ / rho_ + alpha_ * (T - T0_) * g_));

    Scalar error = uEqn_.solve();

    for (const Cell &cell: grid_->localCells())
        u_(cell) += timeStep / rho_ * gradP_(cell);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepBoussinesq::solveTEqn(Scalar timeStep)
{
    T.savePreviousTimeStep(timeStep, 1);

    TEqn_ = (fv::ddt(T, timeStep) + fv::div(u_, T, 0.5)
             == fv::laplacian(kappa_, T, 0.5));

    Scalar error = TEqn_.solve();

    grid_->sendMessages(T);
    T.interpolateFaces();

    return error;
}
