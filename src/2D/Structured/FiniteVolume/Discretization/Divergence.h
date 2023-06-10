#ifndef PHASE_DIVERGENCE_H
#define PHASE_DIVERGENCE_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"
#include "FiniteVolume/Field/VectorField.h"

namespace fv {

template <class T>
FiniteVolumeEquation<T> div(const VectorField &u, Field<T> &phi) {
  FiniteVolumeEquation<T> eqn(phi);

  for (const Cell &cell : phi.grid()->localCells()) {
    for (const InteriorFaceStencil &st : cell.interiorFaceStencils()) {
    }
  }

  return eqn;
}

} // namespace fv

#endif
