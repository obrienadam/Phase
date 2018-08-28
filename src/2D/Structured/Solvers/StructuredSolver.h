#ifndef PHASE_SOLVER_H
#define PHASE_SOLVER_H

#include "System/Input.h"
#include "System/SolverInterface.h"

#include "Structured/StructuredGrid2D/StructuredGrid2D.h"

class Solver: public SolverInterface
{
public:

    Solver(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid);

    virtual Scalar solve(Scalar timeStep) = 0;

protected:

    std::shared_ptr<StructuredGrid2D> grid_;
};


#endif
