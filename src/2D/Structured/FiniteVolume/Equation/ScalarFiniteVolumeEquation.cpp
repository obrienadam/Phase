#include "FiniteVolume/Field/ScalarField.h"

#include "FiniteVolumeEquation.h"

template<>
FiniteVolumeEquation<Scalar>::FiniteVolumeEquation(ScalarField &field)
    :
      Equation(field.grid()->localCells().size(), 10),
      _field(field)
{

}
