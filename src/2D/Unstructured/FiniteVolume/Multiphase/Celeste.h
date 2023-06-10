#ifndef PHASE_CELESTE_H
#define PHASE_CELESTE_H

#include "Math/Matrix.h"
#include "System/Input.h"

#include "SurfaceTensionForce.h"

class Celeste : public SurfaceTensionForce {
public:
  Celeste(const Input &input,
          const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
          const std::shared_ptr<CellGroup> &fluidCells);

  virtual void computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma,
                                          const ScalarGradient &gradGamma);

  virtual void computeInterfaceForces(const ScalarFiniteVolumeField &gamma,
                                      const ScalarGradient &gradGamma);

protected:
  class Stencil {
  public:
    Stencil() {}

    Stencil(const Cell &cell, bool weighted = false);

    void init(bool weighted = false);

    void reset();

    const Cell &cell() const { return *cellPtr_; }

    bool weighted() const { return weighted_; }

    virtual Vector2D grad(const ScalarFiniteVolumeField &phi) const;

    virtual Scalar div(const VectorFiniteVolumeField &u) const;

    virtual Scalar axiDiv(const VectorFiniteVolumeField &u) const;

    virtual Scalar kappa(const VectorFiniteVolumeField &n) const;

  protected:
    virtual void initMatrix();

    static Matrix b_;

    const Cell *cellPtr_ = nullptr;

    bool weighted_;

    Matrix pInv_;

    std::vector<Ref<const Cell>> cells_;

    std::vector<Ref<const Face>> faces_;
  };

  void computeGradGammaTilde(const ScalarFiniteVolumeField &gamma);

  virtual void computeCurvature();

  virtual void computeStencils();

  std::vector<Stencil> kappaStencils_, gradGammaTildeStencils_;
};

#endif
