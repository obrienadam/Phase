#ifndef PHASE_GHOST_CELL_IMMERSED_BOUNDARY_BC_STENCIL_H
#define PHASE_GHOST_CELL_IMMERSED_BOUNDARY_BC_STENCIL_H

#include "Geometry/Tensor2D.h"
#include "GhostCellImmersedBoundary.h"
#include "Math/StaticMatrix.h"

class GhostCellImmersedBoundary::BcStencil {
public:
  BcStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj,
            const FiniteVolumeGrid2D &grid);

  BcStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj,
            const FiniteVolumeGrid2D &grid, const Vector2D &r);

  const Cell &cell() const { return cell_; }

  const std::vector<Ref<const Cell>> &cells() const { return cells_; }

  const std::vector<Scalar> &coeffs() const { return coeffs_; }

  const Point2D &ip() const { return ip_; }

  const Point2D &bp() const { return bp_; }

  const Vector2D &nw() const { return nw_; }

  Scalar length() const { return (ip_ - cell_.get().centroid()).mag(); }

  Scalar bpValue(const ScalarFiniteVolumeField &phi) const;

  Scalar ipValue(const ScalarFiniteVolumeField &phi) const;

  Tensor2D bpGrad(const VectorFiniteVolumeField &u) const;

protected:
  void initCoeffMatrix(const FiniteVolumeGrid2D &grid);

  Ref<const Cell> cell_;

  std::vector<Ref<const Cell>> cells_;

  std::vector<Scalar> coeffs_;

  Point2D ip_, bp_, nw_;

  StaticMatrix<4, 4> Ainv_;
};

#endif
