#ifndef PHASE_SECOND_ORDER_BACKWARDS_TIME_DERIVATIVE_H
#define PHASE_SECOND_ORDER_BACKWARDS_TIME_DERIVATIVE_H

#include "Math/TaylorSeries.h"

#include "FiniteVolumeEquation.h"

namespace fv {
template <typename T>
FiniteVolumeEquation<T> ddt2(FiniteVolumeField<T> &field, Scalar timeStep,
                             const CellGroup &cells) {
  FiniteVolumeEquation<T> eqn(field);

  const auto &field0 = field.oldField(0);
  const auto &field1 = field.oldField(1);

  TaylorSeries ts({0., -field.oldTimeStep(0), -field.oldTimeStep(1)}, 2);

  for (const Cell &cell : cells) {
    eqn.add(cell, cell, ts.coeffs()[0] * cell.volume());
    eqn.addSource(cell, ts.coeffs()[1] * cell.volume() * field0(cell));
    eqn.addSource(cell, ts.coeffs()[2] * cell.volume() * field1(cell));
  }

  return eqn;
}

template <typename T>
FiniteVolumeEquation<T> ddt2(FiniteVolumeField<T> &field, Scalar timeStep) {
  return ddt2(field, timeStep, field.cells());
}
} // namespace fv

#endif
