#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, ScalarFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(grid.nActiveCells(), grid.nActiveCells(), 5),
      b_(grid.nActiveCells()),
      field_(field)
{

}

template<>
Equation<VectorFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, VectorFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(2*grid.nActiveCells(), 2*grid.nActiveCells(), 5),
      b_(2*grid.nActiveCells()),
      field_(field)
{

}
