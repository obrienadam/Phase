#include "Laplacian.h"

namespace fv {

template <>
FiniteVolumeEquation<Vector2D>
laplacian(Scalar gamma, VectorFiniteVolumeField &phi, Scalar theta) {
  FiniteVolumeEquation<Vector2D> eqn(phi);
  const VectorFiniteVolumeField &phi0 = phi.oldField(0);

  for (const Cell &cell : phi.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Scalar coeff =
          gamma * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
      eqn.add(cell, nb.cell(), theta * coeff);
      eqn.add(cell, cell, theta * -coeff);
      eqn.addSource(cell,
                    (1. - theta) * coeff * (phi0(nb.cell()) - phi0(cell)));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Scalar coeff =
          gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

      switch (phi.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED:
        eqn.add(cell, cell, theta * -coeff);
        eqn.addSource(cell, theta * coeff * phi(bd.face()));
        eqn.addSource(cell,
                      (1. - theta) * coeff * (phi0(bd.face()) - phi0(cell)));
        break;

      case VectorFiniteVolumeField::NORMAL_GRADIENT:
        break;
      case VectorFiniteVolumeField::SYMMETRY: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

        eqn.add(cell, cell, theta * -coeff);
        eqn.add(cell, cell, theta * coeff * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * coeff *
                                (dot(phi0(cell), tw) * tw - phi0(cell)));
      } break;

      case VectorFiniteVolumeField::PARTIAL_SLIP: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();
        Scalar lambda = phi.boundaryRefValue(bd.face()).x;

        Scalar a = lambda != 0. ? lambda * coeff / (lambda * coeff - 1.) : 0.;

        eqn.add(cell, cell, theta * -coeff);
        eqn.add(cell, cell, theta * a * coeff * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * coeff *
                                (a * dot(phi0(cell), tw) * tw - phi0(cell)));
      }

      default:
        throw Exception("fv", "laplacian<Vector2D>",
                        "unrecognized or unspecified boundary type.");
      }
    }
  }

  return eqn;
}

template <>
FiniteVolumeEquation<Vector2D> laplacian(const ScalarFiniteVolumeField &gamma,
                                         VectorFiniteVolumeField &phi,
                                         Scalar theta) {
  FiniteVolumeEquation<Vector2D> eqn(phi);
  const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);
  const VectorFiniteVolumeField &phi0 = phi.oldField(0);

  for (const Cell &cell : phi.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Scalar coeff = gamma(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) /
                     nb.rCellVec().magSqr();
      Scalar coeff0 = gamma0(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) /
                      nb.rCellVec().magSqr();
      eqn.add(cell, cell, theta * -coeff);
      eqn.add(cell, nb.cell(), theta * coeff);
      eqn.addSource(cell,
                    (1. - theta) * coeff0 * (phi0(nb.cell()) - phi0(cell)));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Scalar coeff = gamma(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) /
                     bd.rFaceVec().magSqr();
      Scalar coeff0 = gamma0(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) /
                      bd.rFaceVec().magSqr();

      switch (phi.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED:
        eqn.add(cell, cell, theta * -coeff);
        eqn.addSource(cell, theta * coeff * phi(bd.face()));
        eqn.addSource(cell,
                      (1. - theta) * coeff0 * (phi0(bd.face()) - phi0(cell)));
        break;

      case VectorFiniteVolumeField::NORMAL_GRADIENT:
        break;
      case VectorFiniteVolumeField::SYMMETRY: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

        eqn.add(cell, cell, theta * -coeff);
        eqn.add(cell, cell, theta * coeff * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * coeff0 *
                                (dot(phi0(cell), tw) * tw - phi0(cell)));
      } break;

      default:
        throw Exception("fv", "laplacian<Vector2D>",
                        "unrecognized or unspecified boundary type.");
      }
    }
  }

  return eqn;
}

} // namespace fv
