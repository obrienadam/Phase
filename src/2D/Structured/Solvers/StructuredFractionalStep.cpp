#include "StructuredFractionalStep.h"

StructuredFractionalStep::StructuredFractionalStep(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid)
        :
        StructuredSolver(input, grid),
        u(grid),
        p(grid)
{

}

Scalar StructuredFractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correct(timeStep);
}

//- Protected methods

void StructuredFractionalStep::solveUEqn(Scalar timeStep)
{

    for(auto j = 1; j < grid_->nNodesJ() - 1; ++j)
        for(auto i = 1; i < grid_->nNodesI() - 1; ++i)
        {
            auto id = grid_->id(i, j);
        }

}

void StructuredFractionalStep::solvePEqn(Scalar timeStep)
{

}

void StructuredFractionalStep::correct(Scalar timeStep)
{

}