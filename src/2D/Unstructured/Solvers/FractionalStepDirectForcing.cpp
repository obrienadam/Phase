#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryObject.h"

#include "FractionalStepDirectForcing.h"
#include "Math/TrilinosAmesosSparseMatrixSolver.h"
#include "Math/Matrix.h"

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
                    //gradP(nb.cell()) = st.gradPExtend(rho_, gradP);
                }
        }
    }

    grid_->sendMessages(u);
}

Scalar FractionalStepDirectForcing::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.) - src::src(gradP / rho_));

    Scalar error = uEqn_.solve();

    fb.fill(Vector2D(0., 0.));

    reconstructVelocity(timeStep);

    //    for (const auto &ibObj: *ib_)
    //    {
    //        for(const Cell& cell: ibObj->ibCells())
    //        {
    //            DirectForcingImmersedBoundaryObject::Stencil st(u, cell, *ibObj);
    //            fb(cell) = (st.uf() - u(cell)) / timeStep;
    //        }

    //        for(const Cell& cell: ibObj->solidCells())
    //            fb(cell) = (ibObj->velocity(cell.centroid()) - u(cell)) / timeStep;
    //    }

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

void FractionalStepDirectForcing::reconstructVelocity(Scalar timeStep)
{
    auto spSolve = std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm());

    FiniteVolumeEquation<Vector2D> eqn(fb);
    eqn.setSparseSolver(spSolve);

    auto ibCells = grid_->globalCellGroup(ib_->ibCells());
    auto solidCells = grid_->globalCellGroup(ib_->solidCells());

    for(const Cell& cell: grid_->localCells())
    {
        if(ibCells.isInGroup(cell))
        {
            std::vector<const Cell*> stCells;
            std::vector<std::pair<Point2D, Vector2D>> compatPts;

            for(const CellLink &nb: cell.cellLinks())
            {
                if(!solidCells.isInGroup(nb.cell()))
                {
                    stCells.push_back(&nb.cell());

                    if(ibCells.isInGroup(nb.cell()))
                    {
                        auto ibObj = ib_->nearestIbObj(nb.cell().centroid());
                        auto bp = ibObj->nearestIntersect(nb.cell().centroid());
                        auto bu = ibObj->velocity(bp);
                        compatPts.push_back(std::make_pair(bp, bu));
                    }
                }
            }

            auto ibObj = ib_->nearestIbObj(cell.centroid());
            auto bp = ibObj->nearestIntersect(cell.centroid());
            auto bu = ibObj->velocity(bp);

            compatPts.push_back(std::make_pair(bp, bu));

            Matrix A(stCells.size() + compatPts.size(), 6), b(stCells.size() + compatPts.size(), 2);

            for(int i = 0; i < stCells.size(); ++i)
            {
                Point2D x = stCells[i]->centroid();
                A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
            }

            for(int i = 0; i < compatPts.size(); ++i)
            {
                Point2D x = compatPts[i].first;
                A.setRow(i + stCells.size(), {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
            }

            Point2D x = cell.centroid();

            Matrix beta = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * pseudoInverse(A);

            eqn.add(cell, cell, -timeStep);

            for(int i = 0; i < stCells.size(); ++i)
            {
                eqn.add(cell, *stCells[i], beta(0, i) * timeStep);
                eqn.addSource(cell, beta(0, i) * u(*stCells[i]));
            }

            for(int i = 0; i < compatPts.size(); ++i)
                eqn.addSource(cell, beta(0, i + stCells.size()) * compatPts[i].second);

            eqn.addSource(cell, -u(cell));
        }
        else if (solidCells.isInGroup(cell))
        {
            eqn.set(cell, cell, -1.);
            eqn.setSource(cell, (ib_->ibObj(cell.centroid())->velocity(cell.centroid()) - u(cell)) / timeStep);
        }
        else
        {
            eqn.set(cell, cell, -1.);
            eqn.setSource(cell, 0.);
        }
    }

    eqn.solve();
}
