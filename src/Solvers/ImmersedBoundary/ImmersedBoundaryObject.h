#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Circle.h"
#include "Equation.h"
#include "ImmersedBoundaryCeleste.h"

class ImmersedBoundaryObject : public Circle
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, CONTACT_ANGLE, PARTIAL_SLIP};

    ImmersedBoundaryObject(const std::string& name,
                           const FiniteVolumeGrid2D& grid,
                           const ImmersedBoundaryCeleste &csf,
                           const Point2D& center = Vector2D(0., 0.), Scalar radius = 0.);

    Point2D boundaryPoint(const Point2D& pt) const;
    Point2D imagePoint(const Point2D& pt) const;

    void setInternalCells();
    void addBoundaryType(const std::string &name, BoundaryType boundaryType);
    void addBoundaryRefValue(const std::string& name, Scalar boundaryRefValue);

    const CellGroup& cells() const { return grid_.cellGroup(name_ + "_cells"); }
    const std::vector< Ref<const Cell> > boundingCells(const Point2D& pt) const;

    BoundaryType boundaryType(const std::string& name) const { return boundaryTypes_.find(name)->second; }
    Scalar boundaryRefValue(const std::string& name) const { return boundaryRefValues_.find(name)->second; }

    const ImmersedBoundaryCeleste& csf() const { return csf_; }

    const std::string& name() const { return name_; }

protected:

    std::string name_;
    const FiniteVolumeGrid2D &grid_;
    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefValues_;

    const ImmersedBoundaryCeleste &csf_;
};

#endif
