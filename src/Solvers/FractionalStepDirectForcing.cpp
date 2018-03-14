#include "FractionalStepDirectForcing.h"
#include "TimeDerivative.h"
#include "Divergence.h"
#include "Laplacian.h"
#include "Source.h"
#include "DirectForcingImmersedBoundaryObject.h"

FractionalStepDirectForcing::FractionalStepDirectForcing(const Input &input)
        :
        FractionalStep(input),
        fb(addVectorField("fb"))
{

}

Scalar FractionalStepDirectForcing::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    fb.fill(Vector2D(0., 0.));


    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.) - src::src(gradP / rho_, fluid_));

    Scalar error = uEqn_.solve();

    for (const auto &ibObj: *ib_)
    {
        for (const Cell &cell: ibObj->ibCells())
        {
            Vector2D uf = DirectForcingImmersedBoundaryObject::Stencil(u, cell, *ibObj).uf();
            fb(cell) = (uf - u(cell)) / timeStep;
        }

        for (const Cell &cell: ibObj->solidCells())
        {
            Vector2D uf = DirectForcingImmersedBoundaryObject::FieldExtensionStencil(u.oldField(0), cell, *ibObj).uf();
            fb(cell) = (uf - u(cell)) / timeStep;
        }
    }

    for(const Cell& cell: fluid_)
        u(cell) = u.oldField(0)(cell);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.) - src::src(gradP / rho_ - fb, fluid_));

    uEqn_.solve();

    for (const Cell &cell: fluid_)
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

    return error;
}