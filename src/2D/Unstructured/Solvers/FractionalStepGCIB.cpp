#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

#include "FractionalStepGCIB.h"

FractionalStepGCIB::FractionalStepGCIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStep(input, grid),
      ib_(input, grid, fluid_)
{
    ib_.updateCells();
}

Scalar FractionalStepGCIB::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    ib_.computeForce(rho_, mu_, u_, p_);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepGCIB::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u_, timeStep) + fv::div(u_, u_, 0.)  + ib_.velocityBcs(u_) == fv::laplacian(mu_ / rho_, u_, 0.5) - src::src(gradP_ / rho_));

    Scalar error = uEqn_.solve();

    for (const Cell &cell: *fluid_)
        u_(cell) += timeStep / rho_ * gradP_(cell);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepGCIB::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p_) + ib_.bcs(p_) == src::div(u_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p_);

    //- Gradient
    p_.setBoundaryFaces();
    gradP_.compute(*fluid_);

    return error;
}
