#ifndef IMMERSED_BOUNDARY_OBJECT_H
#define IMMERSED_BOUNDARY_OBJECT_H

#include "Shape2D.h"
#include "Equation.h"

class ImmersedBoundaryObject
{
public:

    enum BoundaryType
    {
        FIXED, NORMAL_GRADIENT
    };

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

    bool isInIb(const Point2D &pt) const
    { return shapePtr_->isInside(pt); }

    Shape2D &shape()
    { return *shapePtr_; }

    const Shape2D &shape() const
    { return *shapePtr_; }

    //- Operations
    std::pair<Point2D, Vector2D> intersectionStencil(const Point2D &ptA,
                                                     const Point2D &ptB) const; // returns a intersection point and the edge normal

    void addBoundaryType(const std::string &name, BoundaryType boundaryType);

    void addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const Vector2D &boundaryRefValue);

    void addBoundaryRefValue(const std::string &name, const std::string &value);

    //- Reference to cell zone for iterating

    const CellZone &cells() const
    { return *cells_; }

    const CellGroup &ibCells() const
    { return *ibCells_; }

    const CellGroup &solidCells() const
    { return *solidCells_; }

    const CellGroup &freshlyClearedCells() const
    { return *freshlyClearedCells_; }

    //- Boundary info
    BoundaryType boundaryType(const std::string &name) const
    { return boundaryTypes_.find(name)->second; }

    template<typename T>
    T getBoundaryRefValue(const std::string &name) const;

    //- Motion info if applicable
    const Vector2D& position() const { return shape().centroid(); }

    virtual Vector2D acceleration() const
    { return Vector2D(0., 0.); }

    virtual Vector2D acceleration(const Point2D &point) const
    { return Vector2D(0., 0.); }

    virtual Vector2D velocity() const
    { return Vector2D(0., 0.); }

    virtual Vector2D velocity(const Point2D &point) const
    { return Vector2D(0., 0.); }

    virtual Vector2D normalForce() const
    { return Vector2D(0., 0.); }

    virtual Vector2D shearForce() const
    { return Vector2D(0., 0.); }

    virtual Vector2D force() const { return normalForce() + shearForce(); }

    virtual void computeNormalForce(const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& p);

    virtual void computeShearForce(const ScalarFiniteVolumeField &mu, const VectorFiniteVolumeField& u);

    //- Update (must to be overriden)
    virtual void update(Scalar timeStep) = 0;

    virtual void updateCells() = 0;

    virtual Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const = 0;

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const = 0;

    const std::string &name() const
    { return name_; }

    Label id() const
    { return id_; }

protected:

    std::string name_;
    Label id_;

    FiniteVolumeGrid2D &grid_;

    CellGroup *ibCells_, *solidCells_, *freshlyClearedCells_;
    CellZone *cells_;

    std::shared_ptr<Shape2D> shapePtr_;

    std::map<std::string, BoundaryType> boundaryTypes_;
    std::map<std::string, Scalar> boundaryRefScalars_;
    std::map<std::string, Vector2D> boundaryRefVectors_;
};

#include "GhostCellImmersedBoundaryObject.h"
#include "ForcingCellImmersedBoundaryObject.h"

#endif
