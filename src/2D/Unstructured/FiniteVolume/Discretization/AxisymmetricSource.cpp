#include "AxisymmetricSource.h"

Vector axi::src::src(const ScalarFiniteVolumeField &phi) {
  Vector phiSrc(phi.grid()->localCells().size());

  for (const Cell &cell : phi.cells())
    phiSrc(phi.indexMap()->local(cell, 0)) = phi(cell) * cell.polarVolume();

  return phiSrc;
}

Vector axi::src::src(const VectorFiniteVolumeField &u) {
  Vector divU(2 * u.grid()->localCells().size());

  for (const Cell &cell : u.cells()) {
    Vector2D tmp = u(cell) * cell.polarVolume();
    divU(u.indexMap()->local(cell, 0)) = tmp.x;
    divU(u.indexMap()->local(cell, 1)) = tmp.y;
  }

  return divU;
}

Vector axi::src::div(const VectorFiniteVolumeField &u) {
  Vector divU(u.grid()->localCells().size());

  for (const Cell &cell : u.cells()) {
    Scalar tmp = 0.;
    for (const InteriorLink &nb : cell.neighbours())
      tmp += dot(u(nb.face()), nb.polarOutwardNorm());

    for (const BoundaryLink &bd : cell.boundaries())
      tmp += dot(u(bd.face()), bd.polarOutwardNorm());

    divU(u.indexMap()->local(cell, 0)) = tmp;
  }

  return divU;
}
