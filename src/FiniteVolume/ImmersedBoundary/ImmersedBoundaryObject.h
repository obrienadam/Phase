#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"
#include "BilinearInterpolation.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, CONTACT_ANGLE, PARTIAL_SLIP};

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string& name,
                           const Point2D& center,
                           Scalar radius,
                           Label id,
                           FiniteVolumeGrid2D& grid);

    ImmersedBoundaryObject(const std::string& name,
                           const std::vector<Point2D>& vertices,
                           Label id,
                           FiniteVolumeGrid2D &grid);

    //- The shape
    Shape2D& shape() { return *shapePtr_; }
    const Shape2D& shape() const { return *shapePtr_; }

    //- Accessing stencil points for an ib cell
    const Point2D& boundaryPoint(const Cell &cell) const { return (stencilPoints_.find(cell.id())->second).first; }
    const Point2D& imagePoint(const Cell &cell) const { return (stencilPoints_.find(cell.id())->second).second; }
    const std::pair<Point2D, Point2D>& stencilPoints(const Cell &cell) const { return stencilPoints_.find(cell.id())->second; }

    //- For interpolating the image point
    const std::vector< Ref<const Cell> >& imagePointCells(const Cell &cell) const { return (imagePointStencils_.find(cell.id())->second).first; }
    const BilinearInterpolation& imagePointInterpolation(const Cell &cell) const { return (imagePointStencils_.find(cell.id())->second).second; }

    //- Interpolate a value to the image point
    Scalar imagePointVal(const Cell &cell, const ScalarFiniteVolumeField& field) const;
    Vector2D imagePointVal(const Cell &ell, const VectorFiniteVolumeField& field) const;

    //- Operations
    std::pair<Point2D, Vector2D> intersectionStencil(const Point2D& ptA, const Point2D& ptB) const; // returns a intersection point and the edge normal

    //- Internal cells and boundaries
    void setInternalCells();
    void addBoundaryType(const std::string &name, BoundaryType boundaryType);
    void addBoundaryRefValue(const std::string& name, Scalar boundaryRefValue);

    //- Reference to cell zone for iterating
    const CellGroup& ibCells() const { return *ibCells_; }
    const CellGroup& solidCells() const { return *solidCells_; }

    //- Boundary info
    BoundaryType boundaryType(const std::string& name) const { return boundaryTypes_.find(name)->second; }
    Scalar boundaryRefValue(const std::string& name) const { return boundaryRefValues_.find(name)->second; }

    const std::string& name() const { return name_; }

    void assignId(Label id) { id_ = id; }
    Label id() const { return id_; }

protected:

    void flagIbCells();
    void constructStencils();

    std::string name_;
    Label id_;

    FiniteVolumeGrid2D &grid_;

    const CellGroup* ibCells_, *solidCells_;
    const CellZone* cells_;
    std::shared_ptr<Shape2D> shapePtr_;

    std::map<Label, std::pair<Vector2D, Vector2D> > stencilPoints_;
    std::map<Label, std::pair< std::vector< Ref<const Cell> >, BilinearInterpolation > > imagePointStencils_;

    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefValues_;
};

#endif
