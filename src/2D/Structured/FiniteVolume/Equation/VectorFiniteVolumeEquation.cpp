#include "FiniteVolume/Field/VectorField.h"

#include "FiniteVolumeEquation.h"

template <>
FiniteVolumeEquation<Vector2D>::FiniteVolumeEquation(VectorField &field,
                                                     int nnz)
    : Equation(2 * field.grid()->localCells().size(), nnz), _field(field) {}

template <>
void FiniteVolumeEquation<Vector2D>::add(const Cell &cell0, const Cell &cell1,
                                         Scalar coeff) {
  addCoeff(2 * cell0.lid(), 2 * cell1.gid(), coeff);
  addCoeff(2 * cell0.lid() + 1, 2 * cell1.gid() + 1, coeff);
}

// template<>
// void FiniteVolumeEquation<Vector2D>::add(const Cell &cell0, const Cell
// &cell1, const Vector2D &coeff)
//{
//     addCoeff(2 * cell0.lid(), 2 * cell1.gid(), coeff.x);
//     addCoeff(2 * cell0.lid() + 1, 2 * cell1.gid() + 1, coeff.y);
// }
