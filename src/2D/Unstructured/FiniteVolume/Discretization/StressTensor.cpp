#include "StressTensor.h"
#include "Geometry/Tensor2D.h"

VectorFvmEquation fv::divSigma(Scalar mu, const ScalarFiniteVolumeField &p,
                               VectorFiniteVolumeField &u, Scalar theta) {
  VectorFvmEquation eqn(u, 5);
  const VectorFiniteVolumeField &u0 = u.oldField(0);

  for (const Cell &cell : u.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Scalar coeff = mu * dot(nb.rc(), nb.sf()) / nb.rc().magSqr();

      eqn.add(cell, nb.cell(), coeff * theta);
      eqn.add(cell, cell, -coeff * theta);
      eqn.addSource(cell,
                    -p(nb.face()) * nb.sf() +
                        coeff * (u0(nb.cell()) - u0(cell)) * (1. - theta));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      eqn.addSource(cell, -p(bd.face()) * bd.sf());

      switch (u.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED: {
        Scalar coeff = mu * dot(bd.rf(), bd.sf()) / bd.rf().magSqr();

        eqn.addSource(cell, coeff * u(bd.face()) * theta);
        eqn.add(cell, cell, -coeff * theta);
        eqn.addSource(cell, coeff * (u0(bd.face()) - u0(cell)) * (1. - theta));
      } break;
      case VectorFiniteVolumeField::NORMAL_GRADIENT:
      case VectorFiniteVolumeField::SYMMETRY:
        break;
      default:
        throw Exception("fv", "divSigma",
                        "unrecognized or unspecified boundary type.");
      }
    }
  }

  return eqn;
}

VectorFvmEquation fv::divSigma(const ScalarFiniteVolumeField &mu,
                               const ScalarFiniteVolumeField &p,
                               VectorFiniteVolumeField &u, Scalar theta) {
  VectorFvmEquation eqn(u, 5);
  const VectorFiniteVolumeField &u0 = u.oldField(0);

  for (const Cell &cell : u.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Scalar coeff = mu(nb.face()) * dot(nb.rc(), nb.sf()) / nb.rc().magSqr();

      eqn.add(cell, nb.cell(), coeff * theta);
      eqn.add(cell, cell, -coeff * theta);

      Tensor2D tau0 = mu(nb.face()) * outer(u0(nb.cell()) - u0(cell),
                                            nb.rc() / nb.rc().magSqr());
      eqn.addSource(cell, -p(nb.face()) * nb.sf() +
                              dot(tau0, nb.sf()) * (1. - theta) +
                              dot(tau0.transpose(), nb.sf()));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      switch (u.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED: {
        Scalar coeff = mu(bd.face()) * dot(bd.rf(), bd.sf()) / bd.rf().magSqr();
        eqn.addSource(cell, coeff * theta);
        eqn.add(cell, cell, -coeff * theta);

        Tensor2D tau0 = mu(bd.face()) * outer(u0(bd.face()) - u0(cell),
                                              bd.rf() / bd.rf().magSqr());
        eqn.addSource(cell, -p(bd.face()) * bd.sf() +
                                dot(tau0, bd.sf()) * (1. - theta) +
                                dot(tau0.transpose(), bd.sf()));
      } break;
      case VectorFiniteVolumeField::NORMAL_GRADIENT:
      case VectorFiniteVolumeField::SYMMETRY:
        eqn.addSource(cell, -p(bd.face()) * bd.sf());
        break;
      default:
        throw Exception("fv", "divSigma",
                        "unrecognized or unspecified boundary type.");
      }
    }
  }

  return eqn;
}
