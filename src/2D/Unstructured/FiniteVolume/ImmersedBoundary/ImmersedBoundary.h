#ifndef PHASE_IMMERSED_BOUNDARY_H
#define PHASE_IMMERSED_BOUNDARY_H

#include "CollisionModel.h"
#include "FiniteVolume/Equation/FiniteVolumeEquation.h"
#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "ImmersedBoundaryObject.h"

class ImmersedBoundary {
public:
  typedef boost::geometry::index::quadratic<8, 4> Parameters;

  struct IndexableGetter {
    typedef boost::geometry::model::box<Point2D> result_type;

    result_type
    operator()(const std::shared_ptr<ImmersedBoundaryObject> &ibObj) const {
      return ibObj->shape().boundingBox();
    }
  };

  struct EqualTo {
    bool operator()(const std::shared_ptr<ImmersedBoundaryObject> &lhs,
                    const std::shared_ptr<ImmersedBoundaryObject> &rhs) const {
      return lhs == rhs;
    }
  };

  enum Type { FLUID_CELLS = 1, IB_CELLS = 2, SOLID_CELLS = 3, FRESH_CELLS = 4 };

  ImmersedBoundary(const Input &input,
                   const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                   const std::shared_ptr<CellGroup> &domainCells);

  //- Domain cells
  void setDomainCells(const std::shared_ptr<CellGroup> &domainCells);

  const std::shared_ptr<CellGroup> &domainCells() const { return domainCells_; }

  //- Grid
  const std::shared_ptr<const FiniteVolumeGrid2D> &grid() const {
    return grid_;
  }

  //- Cell zones
  CellGroup ibCells() const;

  CellGroup solidCells() const;

  //- Immersed boundary object access
  virtual std::shared_ptr<ImmersedBoundaryObject> ibObj(const Point2D &pt);

  virtual std::shared_ptr<const ImmersedBoundaryObject>
  ibObj(const Point2D &pt) const;

  virtual std::shared_ptr<ImmersedBoundaryObject> ibObj(const Cell &cell);

  virtual std::shared_ptr<const ImmersedBoundaryObject>
  ibObj(const Cell &cell) const;

  virtual std::shared_ptr<const ImmersedBoundaryObject>
  nearestIbObjSurface(const Cell &cell) const;

  const std::vector<std::shared_ptr<const ImmersedBoundaryObject>> &
  findAllIbObjs(const Point2D &pt) const;

  const std::vector<std::shared_ptr<const ImmersedBoundaryObject>> &
  findAllIbObjs(const Circle &c) const;

  std::shared_ptr<const ImmersedBoundaryObject>
  nearestIbObj(const Point2D &pt) const;

  std::pair<std::shared_ptr<const ImmersedBoundaryObject>, Point2D>
  nearestIntersect(const Point2D &pt) const;

  std::shared_ptr<const ImmersedBoundaryObject>
  ibObj(const std::string &name) const;

  const std::vector<std::shared_ptr<ImmersedBoundaryObject>> &ibObjs() const {
    return ibObjs_;
  }

  std::vector<std::shared_ptr<ImmersedBoundaryObject>>::const_iterator
  begin() const {
    return ibObjs_.begin();
  }

  std::vector<std::shared_ptr<ImmersedBoundaryObject>>::const_iterator
  end() const {
    return ibObjs_.end();
  }

  //- Updates
  virtual void updateIbPositions(Scalar timeStep);

  virtual void updateCells() = 0;

  //- Boundary conditions
  template <class T>
  void copyBoundaryConditions(const FiniteVolumeField<T> &srcField,
                              const FiniteVolumeField<T> &destField) {
    for (auto &ibObj : ibObjs_)
      ibObj->addBoundaryCondition(destField.name(),
                                  ibObj->bcType(srcField.name()),
                                  ibObj->bcRefValue<T>(srcField.name()));
  }

  FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

  void clearFreshCells();

  bool isIbCell(const Cell &cell) const;

  virtual void applyHydrodynamicForce(Scalar rho, Scalar mu,
                                      const VectorFiniteVolumeField &u,
                                      const ScalarFiniteVolumeField &p,
                                      const Vector2D &g = Vector2D(0., 0.));

  virtual void applyHydrodynamicForce(const ScalarFiniteVolumeField &rho,
                                      const ScalarFiniteVolumeField &mu,
                                      const VectorFiniteVolumeField &u,
                                      const ScalarFiniteVolumeField &p,
                                      const Vector2D &g = Vector2D(0., 0.));

  virtual void applyCollisionForce(bool add = false);

  const std::shared_ptr<FiniteVolumeField<int>> &cellStatus() {
    return cellStatus_;
  }

protected:
  static std::shared_ptr<ImmersedBoundaryObject>
  createIbObj(const std::string &name,
              const std::unordered_map<std::string, std::string> &properties);

  mutable std::vector<std::shared_ptr<const ImmersedBoundaryObject>> query_;

  void setCellStatus();

  std::shared_ptr<CellGroup> domainCells_;

  std::shared_ptr<FiniteVolumeField<int>> cellStatus_;

  std::shared_ptr<const FiniteVolumeGrid2D> grid_;

  std::vector<std::shared_ptr<ImmersedBoundaryObject>> ibObjs_;

  //- Fast searching
  boost::geometry::index::rtree<std::shared_ptr<ImmersedBoundaryObject>,
                                Parameters, IndexableGetter, EqualTo>
      rTree_;

  //- Collision model
  std::shared_ptr<CollisionModel> collisionModel_;
};

#endif
