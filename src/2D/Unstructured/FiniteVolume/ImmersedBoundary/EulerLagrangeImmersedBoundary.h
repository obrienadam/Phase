#ifndef PHASE_EULER_LAGRANGE_IMMERSED_BOUNDARY_H
#define PHASE_EULER_LAGRANGE_IMMERSED_BOUNDARY_H

#include "System/NotImplementedException.h"

#include "ImmersedBoundary.h"

class EulerLagrangeImmersedBoundary : public ImmersedBoundary {
public:
  EulerLagrangeImmersedBoundary(
      const std::string &name,
      const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
      const std::shared_ptr<CellGroup> &domainCells);

  void initCircle(const Point2D &center, Scalar radius);

  void updateCells();

  //- These methods aren't used, this is strictly a corrective method
  FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &field) const {
    throw NotImplementedException("EulerLagrangeImmersedBoundaryObject", "bcs");
  }

  FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &field) const {
    throw NotImplementedException("EulerLagrangeImmersedBoundaryObject", "bcs");
  }

  FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

  void computeForce(Scalar rho, Scalar mu, const VectorFiniteVolumeField &u,
                    const ScalarFiniteVolumeField &p,
                    const Vector2D &g = Vector2D(0., 0.)) {}

  void correctVelocity(VectorFiniteVolumeField &u) const;

  Scalar kernel(const Point2D &x, const Point2D &xl) const;

private:
  void initLagrangePoints(int nLagrangePoints);

  void updateLagrangePoints() { initLagrangePoints(lagrangePoints_.size()); }

  // std::shared_ptr<TrilinosAmesosSparseMatrixSolver> solver_;

  Scalar h_;

  std::vector<Point2D> lagrangePoints_;

  std::vector<std::vector<Ref<const Cell>>> lagrangeStencils_;
};

#endif
