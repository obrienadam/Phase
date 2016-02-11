#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, ScalarFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(grid.nCells(), grid.nCells(), 5),
      field_(field)
{

}

template<>
Equation<VectorFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, VectorFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(2*grid.nCells(), 2*grid.nCells(), 5),
      field_(field)
{

}
