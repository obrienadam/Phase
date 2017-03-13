#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"
#include "Interpolation.h"
#include "GhostCellStencil.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT, CONTACT_ANGLE, PARTIAL_SLIP
    };
    enum FlaggingCriterion
    {
        CELL_CENTER_WITHIN, CELL_OVERLAPPING
    };

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string &name,
                           Label id,
                           FiniteVolumeGrid2D &grid);

    //- The shape
    void initCircle(const Point2D &center, Scalar radius);

    void initPolygon(const std::vector<Point2D> &vertices);

    Shape2D &shape()
    { return *shapePtr_; }

    const Shape2D &shape() const
    { return *shapePtr_; }

    //- Stencils
    const std::vector<GhostCellStencil> &stencils() const
    { return stencils_; }

    //- Operations
    std::pair<Point2D, Vector2D> intersectionStencil(const Point2D &ptA,
                                                     const Point2D &ptB) const; // returns a intersection point and the edge normal

    //- Internal cells and boundaries
    void setInternalCells();

    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

    void addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const Vector2D &boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const std::string &value);

    //- Reference to cell zone for iterating
    const CellGroup &ibCells() const
    { return *ibCells_; }

    const CellGroup &solidCells() const
    { return *solidCells_; }

    const CellGroup &freshlyClearedCells() const
    { return *freshlyClearedCells_; }

    const CellGroup &cutCells() const
    { return *cutCells_; }

    //- Boundary info
    BoundaryType boundaryType(const std::string &name) const
    { return boundaryTypes_.find(name)->second; }

    template<typename T>
    T getBoundaryRefValue(const std::string &name) const;

    //- Motion info
    virtual Vector2D velocity(const Point2D &point) const
    { return Vector2D(0., 0.); }

    virtual Vector2D trajectory() const
    { return Vector2D(0., 0.); }

    virtual Scalar angularVelocity() const
    { return 0.; }

    //- Update (meant to be overriden
    virtual void update(Scalar timeStep)
    {} //- By default, do nothing
    virtual void updateCells();

    const std::string &name() const
    { return name_; }

    void assignId(Label id)
    { id_ = id; }

    Label id() const
    { return id_; }

protected:

    void constructStencils();

    std::string name_;
    Label id_;

    FiniteVolumeGrid2D &grid_;

    CellGroup *ibCells_, *solidCells_, *freshlyClearedCells_, *cutCells_;
    CellZone *cells_;
    std::shared_ptr<Shape2D> shapePtr_;

    std::vector<GhostCellStencil> stencils_;

    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefScalars_;
    std::map<std::string, Vector2D> boundaryRefVectors_;
};

#endif
