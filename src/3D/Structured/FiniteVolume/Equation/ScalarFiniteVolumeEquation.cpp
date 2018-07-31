#include "ScalarFiniteVolumeEquation.h"

template<>
ScalarFiniteVolumeEquation::FiniteVolumeEquation(Field<Scalar> &field)
    :
      Equation(field.grid()->nCells(), 15),
      _field(field)
{

}

template<>
void ScalarFiniteVolumeEquation::add(const Cell &cell, const Cell &nb, Scalar coeff)
{
    Equation::addCoeff(cell.id(), nb.id(), coeff);
}

template<>
void ScalarFiniteVolumeEquation::addRhs(const Cell &cell, Scalar val)
{
    Equation::addRhs(cell.id(), val);
}
