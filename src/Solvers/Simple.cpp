#include "Simple.h"
#include "Exception.h"

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