#ifndef PHASE_CELESTE_IMMERSED_BOUNDARY_H
#define PHASE_CELESTE_IMMERSED_BOUNDARY_H

#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"
#include "Geometry/StaticPolyLine2D.h"

#include "Celeste.h"

class CelesteImmersedBoundary : public Celeste {
public:
  class ContactLineStencil {
  public:
    ContactLineStencil(const ImmersedBoundaryObject &ibObj, const Point2D &pt,
                       Scalar theta, const ScalarFiniteVolumeField &gamma);

    void init();

    bool isValid() const {
      return ibObj_ && ((cellA_ && cellB_) || (cellA_ && face_));
    }

    const Cell &cellA() const { return *cellA_; }

    const Cell &cellB() const { return *cellB_; }

    const Face &face() const { return *face_; }

    Scalar theta() const { return gamma_; }

    Scalar alpha() const { return gamma_; }

    Scalar gamma() const { return gamma_; }

    const StaticPolyLine2D<3> &cl() const { return cl_; }

    const Vector2D &ncl() const { return ncl_; }

    Vector2D tcl() const { return (cl_[2] - cl_[0]).unitVec(); }

    template <class T> T interpolate(const FiniteVolumeField<T> &field) const {
      return cellB_ ? alpha_ * field(*cellA_) + (1. - alpha_) * field(*cellB_)
                    : alpha_ * field(*cellA_) + (1. - alpha_) * field(*face_);
    }

  protected:
    static std::queue<Ref<const Cell>> cellQueue_;

    static std::unordered_set<Label> cellIdSet_;

    void init(const ScalarFiniteVolumeField &gamma);

    void init(const Ray2D &r1, const Ray2D &r2,
              const ScalarFiniteVolumeField &gamma);

    void findStencilCells(const Ray2D &r);

    const ImmersedBoundaryObject *ibObj_;

    const Cell *cellA_, *cellB_;

    const Face *face_;

    Scalar theta_, alpha_, gamma_;

    StaticPolyLine2D<3> cl_;

    Vector2D ns_, ncl_;
  };

  CelesteImmersedBoundary(const Input &input,
                          const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                          const std::shared_ptr<CellGroup> &fluidCells,
                          const std::weak_ptr<const ImmersedBoundary> &ib);

  virtual void
  computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma,
                             const ScalarGradient &gradGamma) override;

  Scalar theta(const ImmersedBoundaryObject &ibObj) const;

  void computeContactLineExtension(ScalarFiniteVolumeField &gamma) const;

  ContactLineStencil
  contactLineStencil(const Point2D &xc,
                     const ScalarFiniteVolumeField &gamma) const;

  void computeContactLineStencils(const ScalarFiniteVolumeField &gamma);

  void applyFluidForces(const ScalarFiniteVolumeField &rho,
                        const ScalarFiniteVolumeField &mu,
                        const VectorFiniteVolumeField &u,
                        const ScalarFiniteVolumeField &p,
                        const ScalarFiniteVolumeField &gamma, const Vector2D &g,
                        DirectForcingImmersedBoundary &ib) const;

  void applyFluidForces(const ScalarFiniteVolumeField &rho,
                        const VectorFiniteVolumeField &fb,
                        const ScalarFiniteVolumeField &gamma, const Vector2D &g,
                        DirectForcingImmersedBoundary &ib) const;

protected:
  void computeInterfaceNormals() override;

  std::weak_ptr<const ImmersedBoundary> ib_;

  std::unordered_map<std::string, Scalar> ibContactAngles_;

  CellGroup contactLineExtensionCells_;

  std::vector<ContactLineStencil> contactLineStencils_;
};

#endif
