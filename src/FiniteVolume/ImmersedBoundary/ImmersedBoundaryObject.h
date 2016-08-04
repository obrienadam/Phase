#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"
#include "BilinearInterpolation.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, CONTACT_ANGLE, PARTIAL_SLIP};

    ImmersedBoundaryObject(const std::string& name,
                           const FiniteVolumeGrid2D& grid,
                           const Point2D& center,
                           Scalar radius);

    ImmersedBoundaryObject(const std::string& name,
                           const FiniteVolumeGrid2D& grid,
                           const std::vector<Point2D>& vertices);

    Shape2D& shape() { return *shapePtr_; }
    const Shape2D& shape() const { return *shapePtr_; }

    const Point2D& boundaryPoint(const Cell &cell) const { return (stencilPoints_.find(cell.id())->second).first; }
    const Point2D& imagePoint(const Cell &cell) const { return (stencilPoints_.find(cell.id())->second).second; }
    const std::pair<Point2D, Point2D>& stencilPoints(const Cell &cell) const { return stencilPoints_.find(cell.id())->second; }

    const std::vector< Ref<const Cell> >& imagePointCells(const Cell &cell) const { return (imagePointStencils_.find(cell.id())->second).first; }
    const BilinearInterpolation& imagePointInterpolation(const Cell &cell) const { return (imagePointStencils_.find(cell.id())->second).second; }

    Scalar imagePointVal(const Cell &cell, const ScalarFiniteVolumeField& field) const;
    Vector2D imagePointVal(const Cell &ell, const VectorFiniteVolumeField& field) const;

    void setInternalCells();
    void addBoundaryType(const std::string &name, BoundaryType boundaryType);
    void addBoundaryRefValue(const std::string& name, Scalar boundaryRefValue);

    const CellGroup& cells() const { return grid_.cellGroup(name_ + "_cells"); }

    BoundaryType boundaryType(const std::string& name) const { return boundaryTypes_.find(name)->second; }
    Scalar boundaryRefValue(const std::string& name) const { return boundaryRefValues_.find(name)->second; }

    const std::string& name() const { return name_; }

protected:

    const std::vector< Ref<const Cell> > boundingCells(const Point2D& pt) const;
    void flagIbCells();
    void constructStencils();

    std::string name_;
    const FiniteVolumeGrid2D &grid_;
    std::shared_ptr<Shape2D> shapePtr_;

    std::map<size_t, std::pair<Vector2D, Vector2D> > stencilPoints_;
    std::map<size_t, std::pair< std::vector< Ref<const Cell> >, BilinearInterpolation > > imagePointStencils_;

    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefValues_;
};

#endif
