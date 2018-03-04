#include "FractionalStepEulerLagrange.h"
#include "Source.h"
#include "EulerLagrangeImmersedBoundaryObject.h"

FractionalStepEulerLagrange::FractionalStepEulerLagrange(const Input &input,
                                                         std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{
    for (const auto &ibObj: *ib_)
    {
        if (ibObj->type() != ImmersedBoundaryObject::EULER_LAGRANGE)
        {
            throw Exception("FractionalStepEulerLagrange",
                            "FractionalStepEulerLagrange",
                            "immersed boundary objects must be of type \"euler-lagrange\".");
        }
    }
}

Scalar FractionalStepEulerLagrange::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0) == fv::laplacian(mu_ / rho_, u, 0.5));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (auto &ibObj: *ib_)
        std::static_pointer_cast<EulerLagrangeImmersedBoundaryObject>(ibObj)->correctVelocity(u);

    grid_->sendMessages(u);

    u.interpolateFaces();

    return error;
}

Scalar FractionalStepEulerLagrange::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p, fluid_) == src::div(u, fluid_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.compute(fluid_);

    return error;
}