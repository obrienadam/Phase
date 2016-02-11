#include "Simple.h"

Simple::Simple(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(grid, "u"),
      gradP(grid, "gradP"),
      p(grid, "p"),
      pCorr(grid, "pCorr"),
      uEqn_(grid, u),
      pCorrEqn_(grid, pCorr)
{

}

Scalar Simple::solve(Scalar timeStep)
{

}
