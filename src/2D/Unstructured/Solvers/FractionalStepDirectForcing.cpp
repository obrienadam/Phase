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
      extEqn_(input, gradP, "extEqn")
{
    for(const auto& ibObj: *ib_)
        if(ibObj->type() != ImmersedBoundaryObject::DIRECT_FORCING)
            throw Exception("FractionalStepDirectForcing",
                            "FractionalStepDirectForcing",
                            "immersed boundary object must be of type \"direct-forcing\".");
}

Scalar FractionalStepDirectForcing::solve(Scalar timeStep)
{
    grid_->comm().printf("Performing field extensions...\n");
    solveExtEqns();

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

Scalar FractionalStepDirectForcing::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.)
             == fv::laplacian(mu_ / rho_, u, 0.) - src::src(gradP / rho_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u); //- velocities on non-local procs may be needed for fb

    reconstructVelocity(timeStep);

    for (const Cell &cell: grid_->cells())
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

void FractionalStepDirectForcing::solveExtEqns()
{
    auto ibCells = grid_->globalCellGroup(ib_->ibCells());
    auto solidCells = grid_->globalCellGroup(ib_->solidCells());

    extEqn_.clear();

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
                        auto ba = ibObj->acceleration(bp);
                        compatPts.push_back(std::make_pair(bp, -rho_ * ba));
                    }
                }
            }

            auto ibObj = ib_->nearestIbObj(cell.centroid());
            auto bp = ibObj->nearestIntersect(cell.centroid());
            auto ba = ibObj->acceleration(bp);

            compatPts.push_back(std::make_pair(bp, -rho_ * ba));

            Matrix A(stCells.size() + compatPts.size(), 6);

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

            extEqn_.add(cell, cell, -1.);

            for(int i = 0; i < stCells.size(); ++i)
                extEqn_.add(cell, *stCells[i], beta(0, i));

            for(int i = 0; i < compatPts.size(); ++i)
                extEqn_.addSource(cell, beta(0, i + stCells.size()) * compatPts[i].second);
        }
        else if (solidCells.isInGroup(cell))
        {
            extEqn_.set(cell, cell, -1.);
            extEqn_.setSource(cell, -rho_ * ib_->ibObj(cell.centroid())->acceleration(cell.centroid()));
        }
        else
        {
            extEqn_.set(cell, cell, -1.);
            extEqn_.setSource(cell, gradP(cell));
        }
    }

    extEqn_.solve();
}

void FractionalStepDirectForcing::reconstructVelocity(Scalar timeStep)
{
    FiniteVolumeEquation<Vector2D> eqn(fb);
    eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm()));

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

            if(stCells.size() + compatPts.size() < 6)
            {


                throw Exception("FractionalStepDirectForcing",
                                "reconstructVelocity",
                                "not enough cells to perform velocity interpolation. Cell id = "
                                + std::to_string(cell.globalId()) + ", proc = " + std::to_string(grid_->comm().rank()));
            }

            Matrix A(stCells.size() + compatPts.size(), 6);

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
