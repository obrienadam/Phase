#include "Math/TrilinosAmesosSparseMatrixSolver.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"
#include "FiniteVolume/Discretization/TimeDerivative.h"
#include "FiniteVolume/Discretization/Divergence.h"
#include "FiniteVolume/Discretization/ExplicitDivergence.h"
#include "FiniteVolume/Discretization/Laplacian.h"
#include "FiniteVolume/Discretization/Source.h"
#include "FiniteVolume/Discretization/StressTensor.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryLeastSquaresQuadraticStencil.h"

#include "FractionalStepDFIB.h"

FractionalStepDFIB::FractionalStepDFIB(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStep(input, grid),
      fb_(*addField<Vector2D>("fb", fluid_)),
      fbEqn_(input, fb_, "fbEqn"),
      //extEqn_(input, gradP_, "extEqn"),
      ib_(std::make_shared<DirectForcingImmersedBoundary>(input, grid, fluid_))
{
    ib_->updateCells();
    addField<int>(ib_->cellStatus());
}

Scalar FractionalStepDFIB::solve(Scalar timeStep)
{
    //    //grid_->comm().printf("Performing field extensions...\n");
    //    //solveExtEqns();

    grid_->comm().printf("Updating IB forces and positions...\n");
    //ib_->applyHydrodynamicForce(rho_, mu_, u_, rho_ * p_, g_);
    ib_->applyHydrodynamicForce(rho_, fb_, g_);
    //ib_->applyCollisionForce(true);
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    //    grid_->comm().printf("Peforming field extension...\n");
    //    solveExtEqns();

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

std::shared_ptr<const ImmersedBoundary> FractionalStepDFIB::ib() const
{
    return ib_;
}

Scalar FractionalStepDFIB::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 2);

    //- explicit predictor
    uEqn_ = (fv::ddt(u_, timeStep) + fv::dive(u_, u_, 0.5)
             == fv::laplacian(mu_ / rho_, u_, 0.) - src::src(gradP_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u_);

    //- semi-implicit corrector
    uEqn_ == fv::laplacian(mu_ / rho_, u_, 0.5) - fv::laplacian(mu_ / rho_, u_, 0.5)
            + ib_->velocityBcs(u_, u_, timeStep);

    u_.savePreviousIteration();
    error = uEqn_.solve();

    fb_.fill(Vector2D(0., 0.), grid_->localCells());
    for(const Cell &c: ib_->localIbCells())
        fb_(c) = (u_(c) - u_.prevIteration()(c)) / timeStep;

    for(const Cell &c: ib_->localSolidCells())
        fb_(c) = (u_(c) - u_.prevIteration()(c)) / timeStep;

    grid_->sendMessages(fb_);

    for(const Cell &c: u_.cells())
        u_(c) += timeStep * gradP_(c);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    return error;
}

void FractionalStepDFIB::solveExtEqns()
{
    FiniteVolumeEquation<Vector2D> eqn(gradP_);
    eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm()));

    for(const Cell &c: grid_->localCells())
    {
        if(ib_->localIbCells().isInSet(c))
        {
            auto st = DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil(c, *ib_);
            auto beta = st.interpolationCoeffs(c.centroid());

            eqn.add(c, c, -1.);

            int i = 0;
            for(const Cell* ptr: st.cells())
                eqn.add(c, *ptr, beta(0, i++));

            for(const auto &cmpt: st.compatPts())
                eqn.addSource(c, -beta(0, i++) * cmpt.acceleration());
        }
        else
        {
            eqn.add(c, c, -1.);
            eqn.addSource(c, gradP_(c));
        }
    }

    eqn.solve();
    grid_->sendMessages(gradP_);
}
