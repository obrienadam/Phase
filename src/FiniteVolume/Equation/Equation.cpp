#include "Equation.h"

Equation::Equation(const FiniteVolumeGrid2D &grid)
    :
      grid_(grid),
      spMat_(grid.nCells(), grid.nCells(), 5)
{

}
