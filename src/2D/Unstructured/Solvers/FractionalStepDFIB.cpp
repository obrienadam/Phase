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
      extEqn_(input, gradP_, "extEqn"),
      ib_(std::make_shared<DirectForcingImmersedBoundary>(input, grid, fluid_))
{
    ib_->updateCells();
    addField<int>(ib_->cellStatus());
}

void FractionalStepDFIB::initialize()
{
    FractionalStep::initialize();
    u_.savePreviousTimeStep(0, 2);
}

Scalar FractionalStepDFIB::solve(Scalar timeStep)
{   
    grid_->comm().printf("Updating IB forces and positions...\n");
    computIbForce(timeStep);
    ib_->applyCollisionForce(true);
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Performing field extensions...\n");
    solveExtEqns();

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
    gradP_.fill(Vector2D(0., 0.), ib_->localSolidCells());
    gradP_.fill(Vector2D(0., 0.), ib_->localIbCells());
    gradP_.sendMessages();

    u_.savePreviousTimeStep(timeStep, 2);
    uEqn_ = (fv::ddt(u_, timeStep) + fv::dive(u_, u_, 0.5)
             == fv::laplacian(mu_ / rho_, u_, 0.) - src::src(gradP_));

    Scalar error = uEqn_.solve();
    u_.sendMessages();

    fbEqn_ = ib_->computeForcingTerm(u_, timeStep, fb_);
    fbEqn_.solve();
    fb_.sendMessages();

    uEqn_ == fv::laplacian(mu_ / rho_, u_, 0.5) - fv::laplacian(mu_ / rho_, u_, 0.) + src::src(fb_);
    uEqn_.solve();

    for(const Cell &c: u_.cells())
        u_(c) += timeStep * gradP_(c);

    u_.sendMessages();
    u_.interpolateFaces();

    return error;
}

void FractionalStepDFIB::solveExtEqns()
{
    extEqn_ = ib_->computeFieldExtension(gradP_);
    extEqn_.solve();
    grid_->sendMessages(gradP_);
}

void FractionalStepDFIB::computIbForce(Scalar timeStep)
{
    for(auto &ibObj: *ib_)
    {
        //- Compute the hydro force from the ib force
        Vector2D fh(0., 0.);
        for(const Cell &c: ibObj->cells())
        {
            fh += (u_(c) - u_.oldField(0)(c)) * c.volume() / timeStep;

            for(const InteriorLink &nb: c.neighbours())
            {
                Scalar flux0 = dot(u_.oldField(0)(nb.face()), nb.outwardNorm()) / 2.;
                Scalar flux1 = dot(u_.oldField(1)(nb.face()), nb.outwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(nb.cell())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(nb.cell());
            }

            for(const BoundaryLink &bd: c.boundaries())
            {
                Scalar flux0 = dot(u_.oldField(0)(bd.face()), bd.outwardNorm()) / 2.;
                Scalar flux1 = dot(u_.oldField(1)(bd.face()), bd.outwardNorm()) / 2.;
                fh += std::max(flux0, 0.) * u_.oldField(0)(c) + std::min(flux0, 0.) * u_.oldField(0)(bd.face())
                        + std::max(flux1, 0.) * u_.oldField(1)(c) + std::min(flux1, 0.) * u_.oldField(1)(bd.face());
            }

            fh -= (fb_(c) + g_) * c.volume();
        }

        fh = grid_->comm().sum(rho_ * fh);

        Vector2D fw = ibObj->rho * ibObj->shape().area() * g_;

        if(grid_->comm().isMainProc())
            std::cout << "Hydrodynamic force = " << fh << "\n"
                      << "Weight = " << fw << "\n"
                      << "Net = " << fh + fw << "\n";


        ibObj->applyForce(fh + fw);
    }
}
