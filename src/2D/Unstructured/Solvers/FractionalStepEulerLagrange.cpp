#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/ImmersedBoundary/EulerLagrangeImmersedBoundaryObject.h"

#include "FractionalStepEulerLagrange.h"

FractionalStepEulerLagrange::FractionalStepEulerLagrange(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{

}

Scalar FractionalStepEulerLagrange::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::divc(u, u, 0.5)
             == fv::laplacian(mu_ / rho_, u, 1.5) + 1. / timeStep * ib_->velocityBcs(u));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

//    for(auto& ibObj: *ib_)
//        std::static_pointer_cast<EulerLagrangeImmersedBoundaryObject>(ibObj)->correctVelocity(u);
//
//    grid_->sendMessages(u);

    u.interpolateFaces();

    return error;
}

Scalar FractionalStepEulerLagrange::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p) == src::div(u));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.compute(fluid_);

    return error;
}