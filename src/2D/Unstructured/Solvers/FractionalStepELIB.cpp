#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/ImmersedBoundary/EulerLagrangeImmersedBoundary.h"

#include "FractionalStepELIB.h"

FractionalStepELIB::FractionalStepELIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{

}

Scalar FractionalStepELIB::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);

//    uEqn_ = (fv::ddt(u, timeStep) + fv::divc(u, u, 0.5)
//             == fv::laplacian(mu_ / rho_, u, 1.5) + 1. / timeStep * ib_->velocityBcs(u));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u_);

//    for(auto& ibObj: *ib_)
//        std::static_pointer_cast<EulerLagrangeImmersedBoundaryObject>(ibObj)->correctVelocity(u);
//
//    grid_->sendMessages(u);

    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepELIB::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p_) == src::div(u_));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p_);

    //- Gradient
    p_.setBoundaryFaces();
    gradP_.compute(*fluid_);

    return error;
}
