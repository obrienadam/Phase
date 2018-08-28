#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/SecondOrderExplicitDivergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

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

    grid_->comm().printf("Updating IB positions and cells...\n");
    ib_->updateIbPositions(timeStep);
    ib_->updateCells();

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    grid_->comm().printf("Computing IB forces...\n");

    ib_->applyHydrodynamicForce(rho_, mu_, u_, p_, g_);
    // ib_->applyHydrodynamicForce(rho_, fb_, g_);

    Vector2D force = Vector2D(0., 0.);
    for(const FaceGroup &patch: grid_->patches())
        for(const Face &f: patch)
        {
            Vector2D sf = f.outwardNorm(f.lCell().centroid());

            if(grid_->localCells().isInSet(f.lCell()))
            {
                //force -= rho_ * dot(u_(f), sf) * sf.unitVec();
                force -= p_(f) * sf;
            }
        }

    force = Vector2D(grid_->comm().sum(force.x), grid_->comm().sum(force.y));

    if(grid_->comm().isMainProc())
        std::cout << "Momentum deficit = " << force << "\n";

    return 0;
}

Scalar FractionalStepDFIB::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 2);

    uEqn_ = (fv::ddt(u_, timeStep) + fv::div2e(u_, u_, 0.5)
             == fv::laplacian(mu_ / rho_, u_, 0.) - src::src(gradP_ / rho_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u_); //- velocities on non-local procs may be needed for fb

    fbEqn_ = ib_->computeForcingTerm(u_, timeStep, fb_);
    fbEqn_.solve();
    grid_->sendMessages(fb_);

    for (const Cell &cell: grid_->cells())
        u_(cell) = u_.oldField(0)(cell);

    uEqn_ = (fv::ddt(u_, timeStep) + fv::div2e(u_, u_, 0.5)
             == fv::laplacian(mu_ / rho_, u_, 0.) - src::src(gradP_ / rho_ - fb_));

    error = uEqn_.solve();

    for(const Cell &c: u_.cells())
        u_(c) += timeStep / rho_ * gradP_(c);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    //    const auto &u0 = u_.oldField(0);
    //    for(const Face &f: grid_->interiorFaces())
    //    {
    //        const Cell &l = f.lCell();
    //        const Cell &r = f.rCell();
    //        Scalar g = f.volumeWeight();

    //        u_(f) = g * (u_(l) - u0(l)) + (1. - g) * (u_(r) - u0(r)) + u0(f);
    //    }

    //    u_.setBoundaryFaces();

    return error;
}

void FractionalStepDFIB::solveExtEqns()
{
    //    //auto ibCells = grid_->globalCellGroup(ib_->ibCells());
    //    auto solidCells = grid_->globalCellGroup(ib_->solidCells());

    //    extEqn_.clear();

    //    for(const Cell& cell: grid_->localCells())
    //    {
    //        if(ibCells.isInGroup(cell))
    //        {
    //            std::vector<const Cell*> stCells;
    //            std::vector<std::pair<Point2D, Vector2D>> compatPts;

    //            for(const CellLink &nb: cell.cellLinks())
    //            {
    //                if(!solidCells.isInGroup(nb.cell()))
    //                {
    //                    stCells.push_back(&nb.cell());

    //                    if(ibCells.isInGroup(nb.cell()))
    //                    {
    //                        auto ibObj = ib_->nearestIbObj(nb.cell().centroid());
    //                        auto bp = ibObj->nearestIntersect(nb.cell().centroid());
    //                        auto ba = ibObj->acceleration(bp);
    //                        compatPts.push_back(std::make_pair(bp, -rho_ * ba));
    //                    }
    //                }
    //            }

    //            auto ibObj = ib_->nearestIbObj(cell.centroid());
    //            auto bp = ibObj->nearestIntersect(cell.centroid());
    //            auto ba = ibObj->acceleration(bp);

    //            compatPts.push_back(std::make_pair(bp, -rho_ * ba));

    //            Matrix A(stCells.size() + compatPts.size(), 6);

    //            for(int i = 0; i < stCells.size(); ++i)
    //            {
    //                Point2D x = stCells[i]->centroid();
    //                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    //            }

    //            for(int i = 0; i < compatPts.size(); ++i)
    //            {
    //                Point2D x = compatPts[i].first;
    //                A.setRow(i + stCells.size(), {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    //            }

    //            Point2D x = cell.centroid();

    //            Matrix beta = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

    //            extEqn_.add(cell, cell, -1.);

    //            for(int i = 0; i < stCells.size(); ++i)
    //                extEqn_.add(cell, *stCells[i], beta(0, i));

    //            for(int i = 0; i < compatPts.size(); ++i)
    //                extEqn_.addSource(cell, beta(0, i + stCells.size()) * compatPts[i].second);
    //        }
    //        else if (solidCells.isInGroup(cell))
    //        {
    //            extEqn_.set(cell, cell, -1.);
    //            extEqn_.setSource(cell, -rho_ * ib_->ibObj(cell.centroid())->acceleration(cell.centroid()));
    //        }
    //        else
    //        {
    //            extEqn_.set(cell, cell, -1.);
    //            extEqn_.setSource(cell, gradP(cell));
    //        }
    //    }

    //    extEqn_.solve();
}
