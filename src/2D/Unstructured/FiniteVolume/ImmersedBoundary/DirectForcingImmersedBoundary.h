#ifndef PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_H
#define PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_H

#include "Geometry/Tensor2D.h"
#include "Math/Matrix.h"
#include "Math/StaticMatrix.h"

#include "ImmersedBoundary.h"

class DirectForcingImmersedBoundary : public ImmersedBoundary {
public:
  class LeastSquaresQuadraticStencil;

  //- Constructors, one for circles, another for polygons
  DirectForcingImmersedBoundary(
      const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
      const std::shared_ptr<CellGroup> &domainCells);

  void updateCells() override;

  FiniteVolumeEquation<Vector2D>
  computeForcingTerm(const VectorFiniteVolumeField &u, Scalar timeStep,
                     VectorFiniteVolumeField &fib) const;

  FiniteVolumeEquation<Vector2D>
  computeForcingTerm(const ScalarFiniteVolumeField &rho,
                     const VectorFiniteVolumeField &u, Scalar timeStep,
                     VectorFiniteVolumeField &fib) const;

  FiniteVolumeEquation<Vector2D>
  computeFieldExtension(VectorFiniteVolumeField &gradP) const;

  FiniteVolumeEquation<Vector2D>
  computeFieldExtension(const ScalarFiniteVolumeField &rho,
                        const VectorFiniteVolumeField &sg,
                        VectorFiniteVolumeField &gradP) const;

  virtual FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &phi) const {
    return FiniteVolumeEquation<Scalar>(phi);
  }

  virtual FiniteVolumeEquation<Vector2D>
  bcs(VectorFiniteVolumeField &phi) const {
    return FiniteVolumeEquation<Vector2D>(phi);
  }

  FiniteVolumeEquation<Vector2D>
  velocityBcs(VectorFiniteVolumeField &u, const VectorFiniteVolumeField &uTilde,
              Scalar timeStep) const;

  FiniteVolumeEquation<Vector2D>
  continuityVelocityBcs(VectorFiniteVolumeField &u,
                        const VectorFiniteVolumeField &uTilde,
                        Scalar timeStep) const;

  FiniteVolumeEquation<Vector2D>
  polarVelocityBcs(VectorFiniteVolumeField &u,
                   const VectorFiniteVolumeField &uTilde,
                   Scalar timeStep) const;

  FiniteVolumeEquation<Vector2D>
  continuityPolarVelocityBcs(VectorFiniteVolumeField &u,
                             const VectorFiniteVolumeField &uTilde,
                             Scalar timeStep) const;

  virtual void
  applyHydrodynamicForce(Scalar rho, Scalar mu,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &p,
                         const Vector2D &g = Vector2D(0., 0.)) override;

  virtual void applyHydrodynamicForce(
      const ScalarFiniteVolumeField &rho, const ScalarFiniteVolumeField &mu,
      const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &p,
      const Vector2D &g = Vector2D(0., 0.)) override;

  void applyHydrodynamicForce(Scalar rho, const VectorFiniteVolumeField &fib,
                              const Vector2D &g = Vector2D(0., 0.));

  void applyHydrodynamicForce(const VectorFiniteVolumeField &fib);

  //- Cell group access

  const CellGroup &localIbCells() const { return localIbCells_; }

  const CellGroup &localSolidCells() const { return localSolidCells_; }

  //    const CellGroup &globalIbCells() const
  //    { return globalIbCells_; }

  //    const CellGroup &globalSolidCells() const
  //    { return globalSolidCells_; }

private:
  CellGroup localIbCells_, localSolidCells_;

  CellGroup globalIbCells_, globalSolidCells_;
};

#endif
