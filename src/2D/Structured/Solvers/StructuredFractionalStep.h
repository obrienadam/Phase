#ifndef PHASE_STRUCTURED_FRACTIONAL_STEP_H
#define PHASE_STRUCTURED_FRACTIONAL_STEP_H

#include "StructuredSolver.h"

class StructuredFractionalStep: public StructuredSolver
{
public:

    StructuredFractionalStep(const Input& input, const std::shared_ptr<StructuredGrid2D> &grid);

    Scalar solve(Scalar timeStep);

};


#endif
