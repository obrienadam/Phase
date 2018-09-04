#include "FiniteVolume/Field/VectorField.h"

#include "FiniteVolumeEquation.h"

template<>
FiniteVolumeEquation<Vector2D>::FiniteVolumeEquation(VectorField &field)
    :
      Equation(2 * field.grid()->localCells().size(), 10),
      _field(field)
{

}
