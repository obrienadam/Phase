#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryObject.h"

#include "FractionalStepDirectForcing.h"

FractionalStepDirectForcing::FractionalStepDirectForcing(const Input &input)
        :
        FractionalStep(input),
        divU(*addScalarField("divU")),
        fb(*addVectorField("fb")),
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
    //grid_->comm().printf("Performing field extension...\n");
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
    pExtEqn_.clear();
    CellGroup gr;

    for (const Cell &cell: fluid_)
    {
        auto ibObj = ib_->ibObj(cell.centroid());
        bool extend = false;

        if (ibObj)
            for (const auto &nb: cell.neighbours())
                if (!ibObj->isInIb(nb.cell()))
                    extend = true;

        if (!extend)
        {
            pExtEqn_.set(cell, cell, 1.);
            pExtEqn_.setSource(cell, -p(cell));
        }
        else
        {
            auto st = DirectForcingImmersedBoundaryObject::PressureFieldExtensionStencil(p, cell, *ibObj);
            pExtEqn_.add(cell, st.cells(), st.coeffs());
            pExtEqn_.setSource(cell, rho_ * dot(ibObj->acceleration(st.bp()), st.n()));

            for(const auto &nb: cell.neighbours())
            {
                if(!ibObj->isInIb(nb.cell()))
                    gr.add(nb.cell());
            }
        }
    }

    pExtEqn_.solve();
    grid_->sendMessages(p);
    gradP.compute(gr);
}

Scalar FractionalStepDirectForcing::solveUEqn(Scalar timeStep)
{
    u.interpolateFaces();
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.5) - src::src(gradP / rho_, fluid_));

    Scalar error = uEqn_.solve();

    fb.fill(Vector2D(0., 0.));
    for (const auto &ibObj: *ib_)
    {
        ibObj->computeBoundaryForcing(u, timeStep, fb);
    }

    for (const Cell &cell: fluid_)
        u(cell) = u.oldField(0)(cell);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.5) - src::src(gradP / rho_ - fb, fluid_));

    uEqn_.solve();

    for (const Cell &cell: fluid_)
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

    return error;
}