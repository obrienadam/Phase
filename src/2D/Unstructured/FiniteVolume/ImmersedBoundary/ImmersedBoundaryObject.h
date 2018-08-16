#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_H

#include <unordered_map>

#include "Geometry/Shape2D.h"
#include "Geometry/Tensor2D.h"
#include "FiniteVolume/Motion/Motion.h"
#include "FiniteVolumeGrid2D/Cell/CellGroup.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryConditionType
    {
        FIXED, NORMAL_GRADIENT, VELOCITY
    };

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string &name);

    //- Body info
    const std::string &name() const
    { return _name; }

    //- Cell methods
    std::vector<Ref<const Cell>> cellsWithin(const CellGroup &domainCells) const
    { return domainCells.itemsWithin(*_shape); }

    template<class const_iterator>
    std::vector<Ref<const Cell>> internalPerimeterCells(const_iterator first, const_iterator last, bool includeDiagonals = false) const;

    template<class const_iterator>
    std::vector<Ref<const Cell>> outerPerimeterCells(const_iterator first, const_iterator last, bool includeDiagonals = false) const;

    //- Cell groups
    const CellGroup &ibCells() const
    { return _ibCells; }

    const CellGroup &solidCells() const
    { return _solidCells; }

    void setIbCells(const CellGroup &ibCells);

    void setSolidCells(const CellGroup &solidCells);

    void clear();

    //- Geometry related methods
    virtual void initCircle(const Point2D &center, Scalar radius);

    virtual void initBox(const Point2D &center, Scalar width, Scalar height);

    template<class const_iterator>
    void initPolygon(const_iterator begin, const_iterator end)
    { _shape = std::unique_ptr<Polygon>(new Polygon(begin, end)); }

    Shape2D &shape()
    { return *_shape; }

    const Shape2D &shape() const
    { return *_shape; }

    //- Tests
    bool isInIb(const Point2D &pt) const
    { return _shape->isInside(pt); }

    bool isInIb(const Cell &cell) const
    { return isInIb(cell.centroid()); }

    //- Operations
    LineSegment2D intersectionLine(const LineSegment2D &ln) const;

    LineSegment2D intersectionLine(const Point2D &ptA, const Point2D &ptB) const;

    Point2D nearestIntersect(const Point2D &pt) const
    { return _shape->nearestIntersect(pt); }

    Point2D nearestIntersect(const Ray2D& ray) const
    { return _shape->intersections(ray)[0]; }

    Vector2D nearestEdgeUnitNormal(const Point2D &pt) const;

    //- Boundary methods
    template<class T>
    void addBoundaryCondition(const std::string &name, BoundaryConditionType bcType, T bcRefValue);

    template<class T>
    void addBoundaryCondition(const std::string &name, const std::string &bcType, T bcRefValue);

    template<class T>
    void addBoundaryCondition(const std::string &name, const std::string &bcType, const std::string &bcRefValue);

    BoundaryConditionType bcType(const std::string &name) const;

    template<class T>
    T bcRefValue(const std::string &name) const;

    //- Reference to cell zone for iterating

    const CellGroup &cells() const
    { return _cells; }

    //- Motion info if applicable

    void setMotion(const std::shared_ptr<Motion> &motion)
    { _motion = motion; }

    std::shared_ptr<Motion> motion()
    { return _motion; }

    bool isMoving() const
    { return (bool) _motion; }

    const Vector2D &position() const
    { return shape().centroid(); }

    Vector2D velocity(const Point2D &point) const;

    Vector2D acceleration(const Point2D &point) const;

    Vector2D velocity() const
    { return velocity(position()); }

    Vector2D acceleration() const
    { return acceleration(position()); }

    Scalar theta() const;

    Scalar omega() const;

    Scalar alpha() const;

    Scalar mass() const
    { return rho * _shape->area(); }

    Scalar momentOfInertia() const
    { return rho * _shape->momentOfInertia(); }

    void applyForce(const Vector2D &force)
    { _force = force; }

    const Vector2D &force() const
    { return _force; }

    Scalar torque() const
    { return _torque; }

    //- Update
    void updatePosition(Scalar timeStep);

    //- Public properties
    Scalar rho = 0.;

protected:

    //- Identification
    std::string _name;

    CellGroup _cells, _ibCells, _solidCells;

    //- Boundary type info
    std::unordered_map<std::string, BoundaryConditionType> _bcTypes;

    //- Boundary ref values
    std::unordered_map<std::string, Scalar> _bcScalarRefValues;

    std::unordered_map<std::string, Vector2D> _bcVectorRefValues;

    //- Geometry
    std::unique_ptr<Shape2D> _shape;

    //- Force info
    Vector2D _force;

    Scalar _torque;

    //- Motion
    std::shared_ptr<Motion> _motion;
};

#include "ImmersedBoundaryObject.tpp"

#endif
