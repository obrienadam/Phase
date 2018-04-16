#include "StructuredFractionalStep.h"

StructuredFractionalStep::StructuredFractionalStep(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid)
        :
        StructuredSolver(input, grid)
{

}

Scalar StructuredFractionalStep::solve(Scalar timeStep)
{

}