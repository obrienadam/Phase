#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace fv {

template <class T> FiniteVolumeEquation<T> lap(Scalar gamma, Field<T> &phi) {
  FiniteVolumeEquation<T> eqn(phi);

  for (const Cell &cell : phi.grid()->localCells()) {
    //        for(const InteriorStencil &st: cell.neighbours())
    //        {

    //        }

    //        for(const BoundaryStencil &st: cell.boundaries())
    //        {

    //        }
  }

  return eqn;
}

} // namespace fv

#endif
