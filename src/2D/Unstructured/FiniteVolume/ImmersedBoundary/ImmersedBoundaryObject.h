#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_H

#include "Geometry/Shape2D.h"
#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"
#include "FiniteVolume/Motion/Motion.h"

class SurfaceTensionForce;

class ImmersedBoundaryObject
{
public:

    enum Type
    {
        GHOST_CELL, STEP, QUADRATIC, HIGH_ORDER, EULER_LAGRANGE, DIRECT_FORCING
    };

    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT, PARTIAL_SLIP
    };

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string &name,
                           const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                            const std::shared_ptr<CellGroup> &solverCells);

    //- Body info
    const std::string &name() const
    { return name_; }

    virtual Type type() const = 0;

    //- grid
    const std::shared_ptr<const FiniteVolumeGrid2D> &grid() const
    { return grid_; }

    //- Solver cells
    const std::shared_ptr<CellGroup> &solverCells() const
    { return solverCells_; }

    void setSolverCells(const std::shared_ptr<CellGroup> &solverCells)
    { solverCells_ = solverCells; }

    //- Geometry related methods
    virtual void initCircle(const Point2D &center, Scalar radius);

    virtual void initBox(const Point2D &center, Scalar width, Scalar height);

    template<class const_iterator>
    void initPolygon(const_iterator begin, const_iterator end)
    { shape_ = std::unique_ptr<Polygon>(new Polygon(begin, end)); }

    Shape2D &shape()
    { return *shape_; }

    const Shape2D &shape() const
    { return *shape_; }

    //-

    //- Tests
    bool isInIb(const Point2D &pt) const
    { return shape_->isInside(pt); }

    template<class T>
    bool isInIb(const T &item) const
    { return shape_->isInside(item.centroid()); }

    //- Set/get primary cell zone

    virtual void clear();

    //- Operations
    LineSegment2D intersectionLine(const LineSegment2D &ln) const;

    LineSegment2D intersectionLine(const Point2D &ptA, const Point2D &ptB) const;

    Point2D nearestIntersect(const Point2D &pt) const
    { return shape_->nearestIntersect(pt); }

    Point2D nearestIntersect(const Ray2D& ray) const
    { return shape_->intersections(ray)[0]; }

    Vector2D nearestEdgeNormal(const Point2D &pt) const;

    //- Boundary methods
    void addBoundary(const std::string &name, BoundaryType bType, Scalar ref);

    void addBoundary(const std::string &name, BoundaryType bType, const Vector2D &ref);

    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

    void addBoundaryType(const std::string &name, const std::string &boundaryType);

    void addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const Vector2D &boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const std::string &value);

    //- Reference to cell zone for iterating

    const CellGroup &cells() const
    { return cells_; }

    const CellGroup &ibCells() const
    { return ibCells_; }

    const CellGroup &solidCells() const
    { return solidCells_; }

    const CellGroup &freshCells() const
    { return freshCells_; }

    //- Boundary info
    BoundaryType boundaryType(const std::string &name) const
    { return boundaryTypes_.find(name)->second; }

    template<typename T>
    T getBoundaryRefValue(const std::string &name) const;

    //- Motion info if applicable

    void setMotion(const std::shared_ptr<Motion> &motion)
    { motion_ = motion; }

    std::shared_ptr<Motion> motion()
    { return motion_; }

    bool isMoving() const
    { return (bool) motion_; }

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

    virtual std::vector<Point2D> forcingPoints() const;

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

    virtual void computeForce(const ScalarFiniteVolumeField &rho,
                              const ScalarFiniteVolumeField &mu,
                              const VectorFiniteVolumeField &u,
                              const ScalarFiniteVolumeField &p,
                              const ScalarFiniteVolumeField &gamma,
                              const SurfaceTensionForce &ft,
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
    virtual FiniteVolumeEquation<Scalar> bcs(ScalarFiniteVolumeField &field) const = 0;

    virtual FiniteVolumeEquation<Vector2D> bcs(VectorFiniteVolumeField &field) const = 0;

    virtual FiniteVolumeEquation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    virtual FiniteVolumeEquation<Scalar> pressureBcs(ScalarFiniteVolumeField &p) const;

    virtual void computeBoundaryForcing(const VectorFiniteVolumeField& u,
                                        Scalar timeStep,
                                        VectorFiniteVolumeField &fb) const;

    void clearFreshCells();

    //- Public properties
    Scalar rho = 0.;

protected:

    //- Identification
    std::string name_;

    //- Grid
    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::shared_ptr<CellGroup> solverCells_;

    CellGroup cells_, ibCells_, solidCells_, freshCells_;

    //- Geometry
    std::unique_ptr<Shape2D> shape_;

    //- Boundary condition info
    std::unordered_map<std::string, BoundaryType> boundaryTypes_;

    std::unordered_map<std::string, Scalar> boundaryRefScalars_;

    std::unordered_map<std::string, Vector2D> boundaryRefVectors_;

    //- Force info
    Vector2D force_;

    Scalar torque_;

    //- Motion
    std::shared_ptr<Motion> motion_;
};

#endif
