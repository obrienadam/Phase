#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "ImmersedBoundaryObject.h"
#include "CollisionModel.h"

class SurfaceTensionForce;

class ImmersedBoundary
{
public:

    enum
    {
        FLUID_CELLS = 1, IB_CELLS = 2, SOLID_CELLS = 3, FRESH_CELLS = 4, DEAD_CELLS = 5, BUFFER_CELLS = 6
    };

    ImmersedBoundary(const Input &input, const std::shared_ptr<FiniteVolumeGrid2D> &grid);

    //- Solver/grid info
    const std::shared_ptr<FiniteVolumeGrid2D>& grid()
    { return grid_; }

    std::shared_ptr<const FiniteVolumeGrid2D> grid() const
    { return grid_; }

    //- Cell zones
    void initCellZones(CellZone &zone);

    const CellZone &zone() const
    { return *zone_; }

    const NodeGroup &fluidNodes() const
    { return fluidNodes_; }

    CellGroup ibCells() const;

    CellGroup solidCells() const;

    //- Immersed boundary object access
    std::shared_ptr<const ImmersedBoundaryObject> ibObj(const Point2D &pt) const;

    std::shared_ptr<const ImmersedBoundaryObject> nearestIbObj(const Point2D& pt) const;

    std::pair<std::shared_ptr<const ImmersedBoundaryObject>, Point2D> nearestIntersect(const Point2D& pt) const;

    const ImmersedBoundaryObject &ibObj(const std::string &name) const;

    const std::vector<std::shared_ptr<ImmersedBoundaryObject>>& ibObjs() const
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

    Equation<Scalar> pressureBcs(Scalar rho, ScalarFiniteVolumeField &p) const;

    Equation<Scalar> contactLineBcs(const SurfaceTensionForce &fst, ScalarFiniteVolumeField &gamma) const;

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

protected:

    void setCellStatus();

    const CellZone *zone_ = nullptr;
    NodeGroup fluidNodes_;


    std::shared_ptr<FiniteVolumeGrid2D> grid_;
    std::vector<std::shared_ptr<ImmersedBoundaryObject>> ibObjs_;

    //- Collision model
    std::shared_ptr<CollisionModel> collisionModel_;
};

#endif
