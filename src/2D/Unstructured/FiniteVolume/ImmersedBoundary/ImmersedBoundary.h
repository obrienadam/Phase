#ifndef PHASE_IMMERSED_BOUNDARY_H
#define PHASE_IMMERSED_BOUNDARY_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Equation/FiniteVolumeEquation.h"
#include "ImmersedBoundaryObject.h"
#include "CollisionModel.h"

class ImmersedBoundary
{
public:

    enum Type
    {
        FLUID_CELLS = 1, IB_CELLS = 2, SOLID_CELLS = 3, FRESH_CELLS = 4
    };

    ImmersedBoundary(const Input &input,
                     const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                     const std::shared_ptr<CellGroup> &domainCells);

    //- Domain cells
    void setDomainCells(const std::shared_ptr<CellGroup> &domainCells);

    const std::shared_ptr<CellGroup> &domainCells() const
    { return domainCells_; }

    //- Grid
    const std::shared_ptr<const FiniteVolumeGrid2D> &grid() const
    { return grid_; }

    //- Cell zones
    CellGroup ibCells() const;

    CellGroup solidCells() const;

    //- Immersed boundary object access
    std::shared_ptr<const ImmersedBoundaryObject> ibObj(const Point2D &pt) const;

    std::shared_ptr<const ImmersedBoundaryObject> nearestIbObj(const Point2D &pt) const;

    std::pair<std::shared_ptr<const ImmersedBoundaryObject>, Point2D> nearestIntersect(const Point2D &pt) const;

    std::shared_ptr<const ImmersedBoundaryObject> ibObj(const std::string &name) const;

    const std::vector<std::shared_ptr<ImmersedBoundaryObject>> &ibObjs() const
    { return ibObjs_; }

    std::vector<std::shared_ptr<ImmersedBoundaryObject>>::const_iterator begin() const
    { return ibObjs_.begin(); }

    std::vector<std::shared_ptr<ImmersedBoundaryObject>>::const_iterator end() const
    { return ibObjs_.end(); }

    //- Updates
    virtual void updateIbPositions(Scalar timeStep);

    virtual void updateCells() = 0;

    //- Boundary conditions
    template<class T>
    void copyBoundaryConditions(const FiniteVolumeField<T> &srcField, const FiniteVolumeField<T> &destField)
    {
        for (auto &ibObj: ibObjs_)
            ibObj->addBoundaryCondition(destField.name(), ibObj->bcType(srcField.name()), ibObj->bcRefValue<T>(srcField.name()));
    }

    FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    void clearFreshCells();

    bool isIbCell(const Cell &cell) const;

    virtual void computeForce(Scalar rho,
                              Scalar mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.));

    virtual void computeForce(const ScalarFiniteVolumeField &rho,
                              const ScalarFiniteVolumeField &mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const Vector2D &g = Vector2D(0., 0.));

    const std::shared_ptr<FiniteVolumeField<int>> &cellStatus()
    { return cellStatus_; }

protected:

    void setCellStatus();

    std::shared_ptr<CellGroup> domainCells_;

    std::shared_ptr<FiniteVolumeField<int>> cellStatus_;

    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::vector<std::shared_ptr<ImmersedBoundaryObject>> ibObjs_;

    //- Collision model
    std::shared_ptr<CollisionModel> collisionModel_;
};

#endif
