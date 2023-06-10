#include "FiniteVolumeEquation.h"

template <>
FiniteVolumeEquation<Scalar>::FiniteVolumeEquation(
    ScalarFiniteVolumeField &field, const std::string &name, int nnz)
    : CrsEquation(field.grid()->localCells().size(), nnz), name(name),
      field_(field) {}

template <>
void FiniteVolumeEquation<Scalar>::set(const Cell &cell, const Cell &nb,
                                       Scalar val) {
  setCoeff(field_.indexMap()->local(cell, 0), field_.indexMap()->global(nb, 0),
           val);
}

template <>
void FiniteVolumeEquation<Scalar>::add(const Cell &cell, const Cell &nb,
                                       Scalar val) {
  addCoeff(field_.indexMap()->local(cell, 0), field_.indexMap()->global(nb, 0),
           val);
}

template <>
void FiniteVolumeEquation<Scalar>::addSource(const Cell &cell, Scalar val) {
  addRhs(field_.indexMap()->local(cell, 0), val);
}

template <>
void FiniteVolumeEquation<Scalar>::scale(const Cell &cell, Scalar val) {
  CrsEquation::scaleRow(field_.indexMap()->local(cell, 0), val);
}

template <>
void FiniteVolumeEquation<Scalar>::setSource(const Cell &cell, Scalar val) {
  setRhs(field_.indexMap()->local(cell, 0), val);
}

template <> void FiniteVolumeEquation<Scalar>::remove(const Cell &cell) {
  Index row = field_.indexMap()->local(cell, 0);
  std::fill(colInd_.begin() + rowPtr_[row], colInd_.begin() + rowPtr_[row + 1],
            -1);
  rhs_(row) = 0.;
}

template <> void FiniteVolumeEquation<Scalar>::relax(Scalar relaxationFactor) {
  throw Exception("FiniteVolumeEquation<Scalar>", "relax", "method removed.");
  //    for (const Cell &cell: field_.grid()->localCells())
  //    {
  //        Index row = field_.indexMap()->local(cell, 0);
  //        Scalar &coeff = coeffRef(row, row);

  //        coeff /= relaxationFactor;
  //        _rhs(row) -= (1. - relaxationFactor) * coeff * field_(cell);
  //    }
}

//- Private

template <> void FiniteVolumeEquation<Scalar>::mapFromSparseSolver() {
  const IndexMap &idxMap = *field_.indexMap();

  for (const Cell &cell : field_.grid()->localCells())
    field_(cell) = solver_->x(idxMap.local(cell, 0));
}

template <> Size FiniteVolumeEquation<Scalar>::getRank() const {
  return field_.grid()->localCells().size();
}
