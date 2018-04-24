#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryObject.h"

#include "FractionalStepDirectForcing.h"

FractionalStepDirectForcing::FractionalStepDirectForcing(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid),
        fb(*addField<Vector2D>("fb")),
        pExtEqn_(input, p, "pExtEqn"),
        uExtEqn_(input, u, "uExtEqn")
{
    for(const auto& ibObj: *ib_)
        if(ibObj->type() != ImmersedBoundaryObject::DIRECT_FORCING)
            throw Exception("FractionalStepDirectForcing",
                            "FractionalStepDirectForcing",
                            "immersed boundary object must be of type \"direct-forcing\".");
}

Scalar FractionalStepDirectForcing::solve(Scalar timeStep)
{
    //grid_->comm().printf("Performing field extensions...\n");
    //solveExtEqns();

    grid_->comm().printf("Updating IB positions...\n");
    ib_->update(timeStep);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    grid_->comm().printf("Computing IB forces...\n");
    ib_->computeForce(rho_, mu_, u, p, g_);

    return 0;
}

void FractionalStepDirectForcing::solveExtEqns()
{
    for(const auto &ibObj: *ib_)
    {
        for(const Cell &cell: ibObj->ibCells())
        {
            for(const CellLink &nb: cell.neighbours())
                if(ibObj->isInIb(nb.cell()))
                {
                    auto st = DirectForcingImmersedBoundaryObject::FieldExtensionStencil(nb.cell(), *ibObj);

                    u(nb.cell()) = st.uExtend(u);
                    p(nb.cell()) = rho_ * st.pExtend(p);
                }
        }
    }

    grid_->sendMessages(u);
    grid_->sendMessages(p);

    gradP.compute(fluid_);
}

Scalar FractionalStepDirectForcing::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.) - src::src(gradP / rho_));

    Scalar error = uEqn_.solve();

    fb.fill(Vector2D(0., 0.));

    for (const auto &ibObj: *ib_)
    {
        for(const Cell& cell: ibObj->ibCells())
        {
            DirectForcingImmersedBoundaryObject::Stencil st(u, cell, *ibObj);
            fb(cell) = (st.uf() - u(cell)) / timeStep;
        }

        for(const Cell& cell: ibObj->solidCells())
            fb(cell) = (ibObj->velocity(cell.centroid()) - u(cell)) / timeStep;
    }

    for (const Cell &cell: fluid_)
        u(cell) = u.oldField(0)(cell);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.5) - src::src(gradP / rho_ - fb));

    error = uEqn_.solve();

    for (const Cell &cell: fluid_)
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

    return error;
}
