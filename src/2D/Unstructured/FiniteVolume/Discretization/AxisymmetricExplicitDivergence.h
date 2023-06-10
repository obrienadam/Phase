#ifndef PHASE_AXISYMMETRIC_EXPLICIT_DIVERGENCE_H
#define PHASE_AXISYMMETRIC_EXPLICIT_DIVERGENCE_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi {
template <class T>
FiniteVolumeEquation<T> dive(const VectorFiniteVolumeField &u,
                             FiniteVolumeField<T> &phi, Scalar theta = 1.) {
  FiniteVolumeEquation<T> eqn(phi, 0);
  const VectorFiniteVolumeField &u0 = u.oldField(0);
  const VectorFiniteVolumeField &u1 = u.oldField(1);

  const FiniteVolumeField<T> &phi0 = phi.oldField(0);
  const FiniteVolumeField<T> &phi1 = phi.oldField(1);

  for (const Cell &cell : phi.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Vector2D sf = nb.polarOutwardNorm();

      Scalar flux0 = dot(u0(nb.face()), sf);
      Scalar flux1 = dot(u1(nb.face()), sf);

      eqn.addSource(cell, (std::max(flux0, 0.) * phi0(cell) +
                           std::min(flux0, 0.) * phi0(nb.cell())) *
                              theta);

      eqn.addSource(cell, (std::max(flux1, 0.) * phi1(cell) +
                           std::min(flux1, 0.) * phi1(nb.cell())) *
                              (1. - theta));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Vector2D sf = bd.polarOutwardNorm();
      Scalar flux0 = dot(u0(bd.face()), sf);
      Scalar flux1 = dot(u1(bd.face()), sf);

      switch (phi.boundaryType(bd.face())) {
      case FiniteVolumeField<T>::FIXED:
        eqn.addSource(cell, flux0 * phi0(bd.face()) * theta);
        eqn.addSource(cell, flux1 * phi1(bd.face()) * (1. - theta));
        break;

      case FiniteVolumeField<T>::NORMAL_GRADIENT:
        eqn.addSource(cell, flux0 * phi0(bd.face()) * theta);
        eqn.addSource(cell, flux1 * phi1(bd.face()) * (1. - theta));
        break;

      case FiniteVolumeField<T>::SYMMETRY:
        break;

      default:
        throw Exception("axi", "dive<T>",
                        "unrecognized or unspecified boundary type.");
      }
    }
  }

  return eqn;
}
} // namespace axi

#endif
