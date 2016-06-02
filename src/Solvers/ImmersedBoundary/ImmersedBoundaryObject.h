#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Circle.h"
#include "Equation.h"
#include "ImmersedBoundaryContinuumSurfaceForce.h"

class ImmersedBoundaryObject : public Circle
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, CONTACT_ANGLE, PARTIAL_SLIP};

    ImmersedBoundaryObject(const FiniteVolumeGrid2D& grid, const ImmersedBoundaryContinuumSurfaceForce &csf);

    Point2D boundaryPoint(const Point2D& pt) const;
    Point2D imagePoint(const Point2D& pt) const;

    void setInternalCells();
    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

    const CellGroup& cells() const { return grid_.cellGroup("ibCells"); }
    const std::vector< Ref<const Cell> > boundingCells(const Point2D& pt) const;

    BoundaryType boundaryType(const std::string& name) const { return boundaryTypes_.find(name)->second; }

    const ImmersedBoundaryContinuumSurfaceForce& csf() const { return csf_; }

protected:

    const FiniteVolumeGrid2D &grid_;
    std::map<std::string, BoundaryType> boundaryTypes_;

    const ImmersedBoundaryContinuumSurfaceForce &csf_;
};

#endif
