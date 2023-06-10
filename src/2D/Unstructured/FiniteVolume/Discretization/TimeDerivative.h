#ifndef PHASE_TIME_DERIVATIVE_H
#define PHASE_TIME_DERIVATIVE_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace fv {
template <typename T>
FiniteVolumeEquation<T> ddt(Scalar rho, FiniteVolumeField<T> &field,
                            Scalar timeStep) {
  FiniteVolumeEquation<T> eqn(field);
  const FiniteVolumeField<T> &field0 = field.oldField(0);

  for (const Cell &cell : field.cells()) {
    eqn.add(cell, cell, rho * cell.volume() / timeStep);
    eqn.addSource(cell, -rho * cell.volume() * field0(cell) / timeStep);
  }

  return eqn;
}

template <typename T>
FiniteVolumeEquation<T> ddt(const ScalarFiniteVolumeField &rho,
                            FiniteVolumeField<T> &field, Scalar timeStep) {
  const ScalarFiniteVolumeField &rho0 = rho.oldField(0);
  const FiniteVolumeField<T> &field0 = field.oldField(0);

  FiniteVolumeEquation<T> eqn(field);

  for (const Cell &cell : field.cells()) {
    eqn.add(cell, cell, rho(cell) * cell.volume() / timeStep);
    eqn.addSource(cell, -rho0(cell) * cell.volume() * field0(cell) / timeStep);
  }

  return eqn;
}

template <typename T>
FiniteVolumeEquation<T> ddt(FiniteVolumeField<T> &field, Scalar timeStep) {
  FiniteVolumeEquation<T> eqn(field);
  const FiniteVolumeField<T> &field0 = field.oldField(0);

  for (const Cell &cell : field.cells()) {
    eqn.add(cell, cell, cell.volume() / timeStep);
    eqn.addSource(cell, -cell.volume() * field0(cell) / timeStep);
  }

  return eqn;
}

template <typename T>
FiniteVolumeEquation<T> ddt(FiniteVolumeField<T> &field, Scalar timeStep,
                            const CellGroup &cells) {
  FiniteVolumeEquation<T> eqn(field);
  const FiniteVolumeField<T> &field0 = field.oldField(0);

  for (const Cell &cell : cells) {
    eqn.add(cell, cell, cell.volume() / timeStep);
    eqn.addSource(cell, -cell.volume() * field0(cell) / timeStep);
  }

  return eqn;
}
} // namespace fv

#endif
