#ifndef SIMPLE_H
#define SIMPLE_H

#include "Solver.h"
#include "Input.h"
#include "FiniteVolumeGrid2D.h"
#include "Equation.h"

class Simple : public Solver
{
public:

    Simple(const FiniteVolumeGrid2D& grid, const Input& input);

private:
};

#endif
