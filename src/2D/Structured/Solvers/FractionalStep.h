#ifndef PHASE_FRACTIONAL_STEP_H
#define PHASE_FRACTIONAL_STEP_H

#include "StructuredSolver.h"

class FractionalStep: public StructuredSolver
{
public:

    Scalar solve(Scalar timeStep);
};


#endif //PHASE_FRACTIONALSTEP_H
