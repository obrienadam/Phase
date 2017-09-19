#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"
#include "Motion.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryType{FIXED, NORMAL_GRADIENT, PARTIAL_SLIP};

    //- Constructors, one for circles, another for polygons
    ImmersedBoundaryObject(const std::string &name,
                           Label id,
                           FiniteVolumeGrid2D &grid);

    //- The shape
    void initCircle(const Point2D &center, Scalar radius);

    void initBox(const Point2D &center, Scalar width, Scalar height);

    template<class const_iterator>
    void initPolygon(const_iterator begin, const_iterator end)
    {
        shapePtr_ = std::shared_ptr<Polygon>(new Polygon(begin, end));
    }

    //- Motion
    void setMotion(std::shared_ptr<Motion> motion);

    std::shared_ptr<Motion> motion()
    { return motion_; }

    bool isInIb(const Point2D &pt) const
    { return shapePtr_->isInside(pt); }

    bool isInIb(const Cell& cell) const
    { return shapePtr_->isInside(cell.centroid()); }

    template <class const_iterator>
    bool noneInIb(const_iterator begin, const_iterator end) const
    {
        for(const_iterator it = begin; it != end; ++it)
            if(isInIb(*it))
                return false;
        return true;
    }

    Shape2D &shape()
    { return *shapePtr_; }

    const Shape2D &shape() const
    { return *shapePtr_; }

    //- Set/get primary cell zone
    void setZone(CellZone& zone);

    void clear();

    CellZone& cellZone()
    { return *fluid_; }

    const CellZone& cellZone() const
    { return *fluid_; }

    //- Operations
    LineSegment2D intersectionLine(const LineSegment2D& ln) const;

    Point2D nearestIntersect(const Point2D& pt) const
    { return shapePtr_->nearestIntersect(pt); }

    Vector2D nearestEdgeNormal(const Point2D& pt) const;

    std::pair<Point2D, Vector2D> intersectionStencil(const Point2D &ptA,
                                                     const Point2D &ptB) const; // returns a intersection point and the edge normal

    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

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
    const Vector2D& position() const { return shape().centroid(); }

    virtual Vector2D acceleration() const;

    virtual Vector2D acceleration(const Point2D &point) const;

    virtual Vector2D velocity() const;

    virtual Vector2D velocity(const Point2D &point) const;

    const Vector2D& force() const
    { return force_; }

    Scalar torque() const
    { return torque_; }

    virtual void computeForce(Scalar rho, Scalar mu, const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& p)
    {
        force_ = Vector2D(0., 0.);
        torque_ = 0.;
    }

    //- Update (must to be overriden)
    virtual void update(Scalar timeStep);

    virtual void updateCells();

    virtual Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const = 0;

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const = 0;

    virtual Equation<Vector2D> solidVelocity(VectorFiniteVolumeField& u) const;

    virtual Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField& gamma, Scalar theta) const
    { return bcs(gamma); }

    void clearFreshCells();

    const std::string &name() const
    { return name_; }

    Label id() const
    { return id_; }

protected:

    std::string name_;
    Label id_;

    FiniteVolumeGrid2D &grid_;

    std::shared_ptr<CellZone::ZoneRegistry> zoneRegistry_; //- Registry for these IB cells only
    CellZone cells_, ibCells_, solidCells_, freshCells_, deadCells_;
    CellZone *fluid_ = nullptr;

    std::shared_ptr<Shape2D> shapePtr_;

    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefScalars_;
    std::map<std::string, Vector2D> boundaryRefVectors_;

    Vector2D force_;
    Scalar torque_;

    std::shared_ptr<Motion> motion_;
};

#endif
