#ifndef PHASE_IMMERSED_BOUNDARY_H
#define PHASE_IMMERSED_BOUNDARY_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

#include "FiniteVolume/Multiphase/SurfaceTensionForce.h"

#include "ImmersedBoundaryObject.h"

#include "CollisionModel.h"

class ImmersedBoundary
{
public:

    enum
    {
        FLUID_CELLS = 1, IB_CELLS = 2, SOLID_CELLS = 3, FRESH_CELLS = 4
    };

    ImmersedBoundary(const Input &input,
                     const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                     const std::shared_ptr<CellGroup> &solverCells);

    //- Solver cells
    const std::shared_ptr<CellGroup> &solverCells() const
    { return solverCells_; }

    void setSolverCells(const std::shared_ptr<CellGroup> &solverCells);

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
    void update(Scalar timeStep);

    //- Boundary conditions
    template<class T>
    void copyBoundaryTypes(const FiniteVolumeField<T> &srcField, const FiniteVolumeField<T> &destField)
    {
        for (auto &ibObj: ibObjs_)
            ibObj->addBoundaryType(destField.name(), ibObj->boundaryType(srcField.name()));
    }

    template<class T>
    Equation<T> bcs(FiniteVolumeField<T> &field) const
    {
        Equation<T> eqn(field);

        for (const auto &ibObj: ibObjs_)
            eqn += ibObj->bcs(field);

        return eqn;
    }

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    Equation<Scalar> pressureBcs(ScalarFiniteVolumeField &p) const;

    void clearFreshCells();

    bool isIbCell(const Cell &cell) const;

    void computeForce(Scalar rho,
                      Scalar mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const ScalarFiniteVolumeField &gamma,
                      const SurfaceTensionForce &ft,
                      const Vector2D &g = Vector2D(0., 0.));

    const std::shared_ptr<FiniteVolumeField<int>> &cellStatus()
    { return cellStatus_; }

protected:

    void setCellStatus();

    std::shared_ptr<CellGroup> solverCells_;

    std::shared_ptr<FiniteVolumeField<int>> cellStatus_;

    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::vector<std::shared_ptr<ImmersedBoundaryObject>> ibObjs_;

    //- Collision model
    std::shared_ptr<CollisionModel> collisionModel_;
};

#endif
