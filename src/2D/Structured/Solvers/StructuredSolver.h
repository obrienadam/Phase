#ifndef PHASE_STRUCTURED_SOLVER_H
#define PHASE_STRUCTURED_SOLVER_H

#include "System/Input.h"
#include "System/SolverInterface.h"

#include "Structured/StructuredGrid2D/StructuredGrid2D.h"

class StructuredSolver: public SolverInterface
{
public:

    StructuredSolver(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep) = 0;

protected:

    std::shared_ptr<StructuredGrid2D> grid_;
};


#endif
