#include "ViscousStressTensor.h"
#include "FiniteVolume/Field/JacobianField.h"
#include "FiniteVolume/Field/ScalarGradient.h"

VectorFvmEquation fv::divTau(const ScalarFiniteVolumeField &mu,
                             VectorFiniteVolumeField &u, Scalar theta) {
  VectorFvmEquation eqn(u);
  const ScalarFiniteVolumeField &mu0 = mu.oldField(0);
  const VectorFiniteVolumeField &u0 = u.oldField(0);

  for (const Cell &cell : u.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Scalar coeff = mu(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) /
                     nb.rCellVec().magSqr();
      Scalar coeff0 = mu0(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) /
                      nb.rCellVec().magSqr();
      eqn.add(cell, cell, theta * -coeff);
      eqn.add(cell, nb.cell(), theta * coeff);
      eqn.addSource(cell, (1. - theta) * coeff0 * (u0(nb.cell()) - u0(cell)));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Scalar coeff = mu(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) /
                     bd.rFaceVec().magSqr();
      Scalar coeff0 = mu0(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) /
                      bd.rFaceVec().magSqr();

      switch (u.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED:
        eqn.add(cell, cell, theta * -coeff);
        eqn.addSource(cell, theta * coeff * u(bd.face()));
        eqn.addSource(cell, (1. - theta) * coeff0 * (u0(bd.face()) - u0(cell)));
        break;

      case VectorFiniteVolumeField::NORMAL_GRADIENT:
        break;
      case VectorFiniteVolumeField::SYMMETRY: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

        eqn.add(cell, cell, theta * -coeff);
        eqn.add(cell, cell, theta * coeff * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * coeff0 *
                                (dot(u0(cell), tw) * tw - u0(cell)));
      } break;

      default:
        throw Exception("fv", "divTau",
                        "unrecognized or unspecified boundary type.");
      }
    }

    //- This part is always handled explicitly
    Tensor2D gradUTrans = JacobianField::computeJacobian(u, cell).transpose();
    Vector2D gradMu = ScalarGradient::computeGradient(mu, cell);
    eqn.addSource(cell, dot(gradUTrans, gradMu) * cell.volume());
  }

  return eqn;
}
