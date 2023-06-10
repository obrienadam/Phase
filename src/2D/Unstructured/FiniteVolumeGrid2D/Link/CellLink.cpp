#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolumeGrid2D/Cell/Cell.h"

#include "CellLink.h"

CellLink::CellLink(const Cell &self, const Cell &other)
    : Link(self), cell_(other) {
  rCellVec_ = cell_.centroid() - self_.centroid();
}

Scalar CellLink::alpha(const Point2D &pt) const {
  Scalar ll = (pt - self_.centroid()).mag();
  Scalar lr = (pt - cell_.centroid()).mag();
  return lr / (ll + lr);
}

Scalar CellLink::linearInterpolate(const FiniteVolumeField<Scalar> &phi,
                                   const Point2D &pt) const {
  Scalar alpha = this->alpha(pt);
  return alpha * phi(self_) + (1. - alpha) * phi(cell_);
}

Vector2D CellLink::linearInterpolate(const FiniteVolumeField<Vector2D> &u,
                                     const Point2D &pt) const {
  Scalar alpha = this->alpha(pt);
  return alpha * u(self_) + (1. - alpha) * u(cell_);
}
