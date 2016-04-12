#include "Piso.h"

Piso::Piso(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Simple(grid, input)
{
    nPCorrections_ = input.caseInput().get<int>("Solver.numPressureCorrections", 1);
}

Scalar Piso::solve(Scalar timeStep)
{
    solveUEqn(timeStep);

    for(int i = 0; i < nPCorrections_; ++i)
    {
        solvePCorrEqn();
        correctPressure();
        correctVelocity();
    }

    return 0;
}
