#include "AxisymmetricStressTensor.h"

VectorFvmEquation axi::divSigma(Scalar mu, const ScalarFiniteVolumeField &p,
                                VectorFiniteVolumeField &u, Scalar theta) {
  FiniteVolumeEquation<Vector2D> eqn(u, 5);
  const VectorFiniteVolumeField &u0 = u.oldField(0);

  for (const Cell &cell : u.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Vector2D sf = nb.polarOutwardNorm();
      Scalar flux = mu * dot(nb.rc(), sf) / nb.rc().magSqr();

      eqn.add(cell, cell, -theta * flux);
      eqn.add(cell, nb.cell(), theta * flux);
      eqn.addSource(cell, -p(nb.face()) * sf +
                              (1. - theta) * flux * (u0(nb.cell()) - u0(cell)));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Vector2D sf = bd.polarOutwardNorm();
      Scalar flux = mu * dot(bd.rf(), sf) / bd.rf().magSqr();

      switch (u.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED:
        eqn.add(cell, cell, -theta * flux);
        eqn.addSource(cell, theta * flux * u0(bd.face()));
        eqn.addSource(cell, (1. - theta) * flux * (u0(bd.face()) - u0(cell)));
        break;

      case VectorFiniteVolumeField::NORMAL_GRADIENT:
        break;
      case VectorFiniteVolumeField::SYMMETRY: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

        eqn.add(cell, cell, -theta * flux);
        eqn.add(cell, cell, theta * flux * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * flux *
                                (dot(u0(cell), tw) * tw - u0(cell)));
      } break;

      default:
        throw Exception("axi", "divSigma",
                        "unrecognized or unspecified boundary type.");
      }

      eqn.addSource(cell, -p(bd.face()) * sf);
    }

    //- Non-polar volume since we just want the rz face
    Scalar r = cell.centroid().x;
    eqn.add(cell, cell, -mu * cell.volume() / r * Vector2D(theta, 0.));
    eqn.addSource(cell, p(cell) * Vector2D(cell.volume(), 0.) -
                            mu * cell.volume() * u0(cell).x / r *
                                Vector2D(1. - theta, 0.));
  }

  return eqn;
}

VectorFvmEquation axi::divSigma(const ScalarFiniteVolumeField &rho,
                                const ScalarFiniteVolumeField &mu,
                                const ScalarFiniteVolumeField &p,
                                VectorFiniteVolumeField &u, Scalar theta) {
  FiniteVolumeEquation<Vector2D> eqn(u, 5);
  const ScalarFiniteVolumeField &mu0 = mu.oldField(0);
  const VectorFiniteVolumeField &u0 = u.oldField(0);

  for (const Cell &cell : u.cells()) {
    for (const InteriorLink &nb : cell.neighbours()) {
      Vector2D sf = nb.polarOutwardNorm();
      Scalar flux = mu(nb.face()) * dot(nb.rc(), sf) / nb.rc().magSqr();
      Scalar flux0 = mu0(nb.face()) * dot(nb.rc(), sf) / nb.rc().magSqr();

      eqn.add(cell, cell, -theta * flux);
      eqn.add(cell, nb.cell(), theta * flux);
      eqn.addSource(cell,
                    -p(nb.face()) * sf * rho(cell) / rho(nb.face()) +
                        (1. - theta) * flux0 * (u0(nb.cell()) - u0(cell)));
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Vector2D sf = bd.polarOutwardNorm();
      Scalar flux = mu(bd.face()) * dot(bd.rf(), sf) / bd.rf().magSqr();
      Scalar flux0 = mu0(bd.face()) * dot(bd.rf(), sf) / bd.rf().magSqr();

      switch (u.boundaryType(bd.face())) {
      case VectorFiniteVolumeField::FIXED:
        eqn.add(cell, cell, -theta * flux);
        eqn.addSource(cell, theta * flux * u0(bd.face()));
        eqn.addSource(cell, (1. - theta) * flux0 * (u0(bd.face()) - u0(cell)));
        break;

      case VectorFiniteVolumeField::NORMAL_GRADIENT:
        break;
      case VectorFiniteVolumeField::SYMMETRY: {
        Vector2D tw = bd.outwardNorm().tangentVec().unitVec();

        eqn.add(cell, cell, -theta * flux);
        eqn.add(cell, cell, theta * flux * outer(tw, tw));
        eqn.addSource(cell, (1. - theta) * flux0 *
                                (dot(u0(cell), tw) * tw - u0(cell)));
      } break;

      default:
        throw Exception("axi", "divSigma",
                        "unrecognized or unspecified boundary type.");
      }

      eqn.addSource(cell, -p(bd.face()) * sf * rho(cell) / rho(bd.face()));
    }

    //- Non-polar volume since we just want the rz face
    Scalar r = cell.centroid().x;
    eqn.add(cell, cell, -mu(cell) * cell.volume() / r * Vector2D(theta, 0.));
    eqn.addSource(cell, p(cell) * Vector2D(cell.volume(), 0.) -
                            mu0(cell) * cell.volume() * u0(cell).x / r *
                                Vector2D(1. - theta, 0.));
  }

  return eqn;
}
