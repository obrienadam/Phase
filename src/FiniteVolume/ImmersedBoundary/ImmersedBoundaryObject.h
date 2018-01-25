#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"
#include "Motion.h"

class ImmersedBoundary;

class ImmersedBoundaryObject
{
public:

    typedef typename std::vector<std::shared_ptr<ImmersedBoundaryObject>>::iterator iterator;
    typedef typename std::vector<std::shared_ptr<ImmersedBoundaryObject>>::const_iterator const_iterator;

    enum Type
    {
        GHOST_CELL, STEP, QUADRATIC, HIGH_ORDER
    };

    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT, PARTIAL_SLIP
    };

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string &name,
                           Label id,
                           const ImmersedBoundary &ib,
                           const std::shared_ptr<FiniteVolumeGrid2D>& grid);

    virtual Type type() const = 0;

    //- The shape
    void initCircle(const Point2D &center, Scalar radius);

    void initBox(const Point2D &center, Scalar width, Scalar height);

    template<class const_iterator>
    void initPolygon(const_iterator begin, const_iterator end)
    {
        shape_ = std::unique_ptr<Polygon>(new Polygon(begin, end));
    }

    Shape2D &shape()
    { return *shape_; }

    const Shape2D &shape() const
    { return *shape_; }

    bool isInIb(const Point2D &pt) const
    { return shape_->isInside(pt); }

    template<class T>
    bool isInIb(const T &item) const
    { return shape_->isInside(item.centroid()); }

    template<class const_iterator>
    bool allInIb(const_iterator begin, const_iterator end) const
    {
        for (const_iterator it = begin; it != end; ++it)
            if (!isInIb(*it))
                return false;
        return true;
    }

    template<class const_iterator>
    bool noneInIb(const_iterator begin, const_iterator end) const
    {
        for (const_iterator it = begin; it != end; ++it)
            if (isInIb(*it))
                return false;
        return true;
    }

    //- Motion
    void setMotion(std::shared_ptr<Motion> motion);

    std::shared_ptr<Motion> motion()
    { return motion_; }

    //- Set/get primary cell zone
    void setZone(CellZone &zone);

    virtual void clear();

    CellZone &cellZone()
    { return *fluid_; }

    const CellZone &cellZone() const
    { return *fluid_; }

    const std::shared_ptr<FiniteVolumeGrid2D> &grid() const
    { return grid_; }

    //- Operations
    LineSegment2D intersectionLine(const LineSegment2D &ln) const;

    LineSegment2D intersectionLine(const Point2D &ptA, const Point2D &ptB) const;

    Point2D nearestIntersect(const Point2D &pt) const
    { return shape_->nearestIntersect(pt); }

    Vector2D nearestEdgeNormal(const Point2D &pt) const;

    std::pair<Point2D, Vector2D> intersectionStencil(const Point2D &ptA,
                                                     const Point2D &ptB) const; // returns a intersection point and the edge normal

    //- Boundary methods
    void addBoundary(const std::string &name, BoundaryType bType, Scalar ref);

    void addBoundary(const std::string &name, BoundaryType bType, const Vector2D& ref);

    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

    void addBoundaryType(const std::string &name, const std::string &boundaryType);

    void addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const Vector2D &boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const std::string &value);

    //- Reference to cell zone for iterating

    const CellZone &cells() const
    { return cells_; }

    const CellZone &ibCells() const
    { return ibCells_; }

    const CellZone &solidCells() const
    { return solidCells_; }

    const CellZone &freshCells() const
    { return freshCells_; }

    const CellZone &deadCells() const
    { return deadCells_; }

    //- Boundary info
    BoundaryType boundaryType(const std::string &name) const
    { return boundaryTypes_.find(name)->second; }

    template<typename T>
    T getBoundaryRefValue(const std::string &name) const;

    //- Motion info if applicable
    bool isMoving() const
    { return (bool) motion_; }

    const Vector2D &position() const
    { return shape().centroid(); }

    Vector2D acceleration() const;

    Vector2D acceleration(const Point2D &point) const;

    Vector2D velocity() const;

    Vector2D velocity(const Point2D &point) const;

    Scalar theta() const;

    Scalar omega() const;

    Scalar alpha() const;

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

    void addForce(const Vector2D &force)
    {
        force_ += force;
    }

    Scalar mass() const
    { return rho * shape_->area(); }

    Scalar momentOfInertia() const
    { return rho * shape_->momentOfInertia(); }

    const Vector2D &force() const
    { return force_; }

    Scalar torque() const
    { return torque_; }

    //- Update
    void update(Scalar timeStep);

    virtual void updateCells();

    //- Boundary conditions
    virtual Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const = 0;

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const = 0;

    virtual Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    virtual Equation<Scalar> pressureBcs(Scalar rho, ScalarFiniteVolumeField &p) const
    { throw Exception("ImmersedBoundaryObject", "pressureBcs", "not implemented."); }

    virtual Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const
    { return bcs(gamma); }

    void clearFreshCells();

    const std::string &name() const
    { return name_; }

    Label id() const
    { return id_; }

    //- Public properties
    Scalar rho = 0.;

protected:

    std::string name_;
    Label id_;

    std::shared_ptr<FiniteVolumeGrid2D> grid_;

    std::shared_ptr<CellZone::ZoneRegistry> zoneRegistry_; //- Registry for these IB cells only
    CellZone cells_, ibCells_, solidCells_, freshCells_, deadCells_;
    CellZone *fluid_ = nullptr;

    const ImmersedBoundary *ib_ = nullptr;

    std::unique_ptr<Shape2D> shape_;

    std::unordered_map<std::string, BoundaryType> boundaryTypes_;
    std::unordered_map<std::string, Scalar> boundaryRefScalars_;
    std::unordered_map<std::string, Vector2D> boundaryRefVectors_;

    Vector2D force_;
    Scalar torque_;

    std::shared_ptr<Motion> motion_;
};

#endif
