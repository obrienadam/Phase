#ifndef SOLVER_H
#define SOLVER_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class Solver
{
public:

    enum TimeDependent{ON, OFF};

    Solver(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual std::string info();

private:

    const FiniteVolumeGrid2D& grid_;

    TimeDependent timeDependent_;
    int maxIterations_;
};

#endif
