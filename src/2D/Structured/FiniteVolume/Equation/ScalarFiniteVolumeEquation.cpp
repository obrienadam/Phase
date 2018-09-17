#include "FiniteVolume/Field/ScalarField.h"

#include "FiniteVolumeEquation.h"

template<>
FiniteVolumeEquation<Scalar>::FiniteVolumeEquation(ScalarField &field, int nnz)
    :
      Equation(field.grid()->localCells().size(), nnz),
      _field(field)
{

}

template<>
void FiniteVolumeEquation<Scalar>::add(const Cell &cell0, const Cell &cell1, Scalar coeff)
{
    addCoeff(cell0.lid(), cell1.gid(), coeff);
}

template<>
void FiniteVolumeEquation<Scalar>::addSource(const Cell &cell, const Scalar &src)
{
    Equation::addRhs(cell.lid(), src);
}
