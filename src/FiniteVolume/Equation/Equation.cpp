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
      field_(field)
{

}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator =(const Term& term)
{
    spMat_.assemble(term.coefficients());
    const auto &sources = term.sources();

    for(int i = 0, end = b_.size(); i < end; ++i)
        b_[i] = sources[i];

    return *this;
}
