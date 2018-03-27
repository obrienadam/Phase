#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

#include "FractionalStepBoussinesq.h"

FractionalStepBoussinesq::FractionalStepBoussinesq(const Input &input)
        :
        FractionalStep(input),
        T(*addScalarField(input, "T")),
        TEqn_(input, T, "TEqn")
{
    alpha_ = input.caseInput().get<Scalar>("Properties.alpha", 0.00369);
    T0_ = input.caseInput().get<Scalar>("Properties.T0", 273.);
    kappa_ = input.caseInput().get<Scalar>("Properties.kappa", 0.02);
}

Scalar FractionalStepBoussinesq::solve(Scalar timeStep)
{
    grid_->comm().printf("Updating IB positions...\n");
    ib_->update(timeStep);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveTEqn(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    grid_->comm().printf("Computing IB forces...\n");
    ib_->computeForce(rho_, mu_, u, p, g_);

    return 0;
}

Scalar FractionalStepBoussinesq::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.5) + ib_->velocityBcs(u)
             == fv::laplacian(mu_ / rho_, u, fluid_, 0.5) - src::src(gradP / rho_ + alpha_ * (T - T0_) * g_, fluid_));

    Scalar error = uEqn_.solve();

    for (const Cell &cell: grid_->localActiveCells())
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

    return error;
}

Scalar FractionalStepBoussinesq::solveTEqn(Scalar timeStep)
{
    T.savePreviousTimeStep(timeStep, 1);

    TEqn_ = (fv::ddt(T, timeStep) + fv::div(u, T, 0.5) + ib_->bcs(T)
             == fv::laplacian(kappa_, T, 0.5));

    Scalar error = TEqn_.solve();

    grid_->sendMessages(T);
    T.interpolateFaces();

    return error;
}