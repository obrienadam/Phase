#include "ScalarGradient.h"

Vector2D ScalarGradient::computeGradient(const ScalarFiniteVolumeField &phi,
                                         const Cell &c) {
  Vector2D sumSf(0., 0.);
  Vector2D tmp(0., 0.);

  for (const InteriorLink &nb : c.neighbours()) {
    Vector2D rc = nb.rCellVec();
    auto gradPhiF = (phi(nb.cell()) - phi(c)) * rc / rc.magSqr();

    Vector2D sf = nb.outwardNorm().abs();
    tmp += pointwise(gradPhiF, sf);
    sumSf += sf;
  }

  for (const BoundaryLink &bd : c.boundaries()) {
    Vector2D rf = bd.rFaceVec();
    auto gradPhiF = (phi(bd.face()) - phi(c)) * rf / rf.magSqr();
    Vector2D sf = bd.outwardNorm().abs();
    tmp += pointwise(gradPhiF, sf);
    sumSf += sf;
  }

  return Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y);
}

ScalarGradient::ScalarGradient(const ScalarFiniteVolumeField &phi,
                               const std::shared_ptr<const CellGroup> &cells)
    : VectorFiniteVolumeField(phi.grid(), "grad" + phi.name(), Vector2D(0., 0.),
                              true, false, cells),
      phi_(phi) {}

void ScalarGradient::computeFaces() {
  VectorFiniteVolumeField &gradPhi = *this;

  for (const Face &face : grid_->interiorFaces()) {
    Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
    gradPhi(face) =
        (phi_(face.rCell()) - phi_(face.lCell())) * rc / rc.magSqr();
  }

  for (const Face &face : grid_->boundaryFaces()) {
    Vector2D rf = face.centroid() - face.lCell().centroid();
    gradPhi(face) = (phi_(face) - phi_(face.lCell())) * rf / rf.magSqr();
  }
}

void ScalarGradient::compute(const CellGroup &group, Method method) {
  computeFaces();
  VectorFiniteVolumeField &gradPhi = *this;

  // std::fill(gradPhi.begin(), gradPhi.end(), Vector2D(0., 0.));

  switch (method) {
  case FACE_TO_CELL:
    for (const Cell &cell : group) {
      Vector2D sum(0., 0.), tmp(0., 0.);

      for (const InteriorLink &nb : cell.neighbours()) {
        Vector2D sf = nb.outwardNorm().abs();
        tmp += pointwise(gradPhi(nb.face()), sf);
        sum += sf;
      }

      for (const BoundaryLink &bd : cell.boundaries()) {
        Vector2D sf = bd.outwardNorm().abs();
        tmp += pointwise(gradPhi(bd.face()), sf);
        sum += sf;
      }

      gradPhi(cell) = Vector2D(tmp.x / sum.x, tmp.y / sum.y);
    }
    break;
  case GREEN_GAUSS_CELL:
    for (const Cell &cell : group) {
      for (const InteriorLink &nb : cell.neighbours()) {
        Scalar g = nb.distanceWeight();
        Scalar phiF = g * phi_(cell) + (1. - g) * phi_(nb.cell());
        gradPhi(cell) += phiF * nb.outwardNorm();
      }

      for (const BoundaryLink &bd : cell.boundaries())
        gradPhi(cell) += phi_(bd.face()) * bd.outwardNorm();

      gradPhi(cell) /= cell.volume();
    }
    break;
  case GREEN_GAUSS_NODE:
    for (const Cell &cell : group) {
      for (const InteriorLink &nb : cell.neighbours()) {
        auto lNodeWeights = nb.face().lNode().distanceWeights();
        auto rNodeWeights = nb.face().rNode().distanceWeights();
        Scalar phiLN = 0, phiRN = 0;
        int i = 0;
        for (const Cell &cell : nb.face().lNode().cells())
          phiLN += lNodeWeights[i++] * phi_(cell);

        i = 0;
        for (const Cell &cell : nb.face().rNode().cells())
          phiRN += rNodeWeights[i++] * phi_(cell);

        Scalar phiF = (phiLN + phiRN) / 2.;
        gradPhi(cell) += phiF * nb.outwardNorm();
      }

      for (const BoundaryLink &bd : cell.boundaries()) {
        auto lNodeWeights = bd.face().lNode().distanceWeights();
        auto rNodeWeights = bd.face().rNode().distanceWeights();
        Scalar phiLN = 0, phiRN = 0;
        int i = 0;
        for (const Cell &cell : bd.face().lNode().cells())
          phiLN += lNodeWeights[i++] * phi_(cell);

        i = 0;
        for (const Cell &cell : bd.face().rNode().cells())
          phiRN += rNodeWeights[i++] * phi_(cell);

        Scalar phiF = (phiLN + phiRN) / 2.;
        gradPhi(cell) += phiF * bd.outwardNorm();
      }

      gradPhi(cell) /= cell.volume();
    }
    break;
  }
}

void ScalarGradient::compute(Method method) { compute(*cellGroup_, method); }

void ScalarGradient::computeAxisymmetric(const CellGroup &cells) {
  computeFaces();
  fill(Vector2D(0., 0.), cells);

  auto &gradPhi = *this;

  for (const Cell &cell : cells) {
    Vector2D sum(0., 0.), tmp(0., 0.);

    for (const InteriorLink &nb : cell.neighbours()) {
      Vector2D sf = nb.polarOutwardNorm().abs();
      tmp += pointwise(gradPhi(nb.face()), sf);
      sum += sf;
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Vector2D sf = bd.polarOutwardNorm().abs();
      tmp += pointwise(gradPhi(bd.face()), sf);
      sum += sf;
    }

    gradPhi(cell) = Vector2D(tmp.x / sum.x, tmp.y / sum.y);
  }
}

void ScalarGradient::computeAxisymmetric(const ScalarFiniteVolumeField &cw,
                                         const ScalarFiniteVolumeField &fw,
                                         const CellGroup &cells) {
  computeFaces();
  fill(Vector2D(0., 0.), cells);

  auto &gradPhi = *this;

  for (const Cell &cell : cells) {
    Vector2D sum(0., 0.), tmp(0., 0.);

    for (const InteriorLink &nb : cell.neighbours()) {
      Vector2D sf = nb.polarOutwardNorm().abs();
      tmp += pointwise(gradPhi(nb.face()) / fw(nb.face()), sf);
      sum += sf;
    }

    for (const BoundaryLink &bd : cell.boundaries()) {
      Vector2D sf = bd.polarOutwardNorm().abs();
      tmp += pointwise(gradPhi(bd.face()) / fw(bd.face()), sf);
      sum += sf;
    }

    gradPhi(cell) = cw(cell) * Vector2D(tmp.x / sum.x, tmp.y / sum.y);
  }
}
