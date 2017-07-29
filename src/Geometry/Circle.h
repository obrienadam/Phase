#ifndef CIRCLE_H
#define CIRCLE_H

#include <math.h>

#include "Types.h"
#include "Shape2D.h"
#include "Polygon.h"

class Circle : public Shape2D
{
public:

    Circle(Point2D center, Scalar radius);

    Type type() const
    { return CIRCLE; }

    void init(const Point2D &center, Scalar radius);

    //- Circle parameters
    const Point2D &centroid() const
    { return center_; }

    Scalar area() const
    { return area_; }

    Scalar momentOfInertia() const
    { return momentOfInertia_; }

    Scalar perimeter() const
    { return 2.*M_PI*radius_; }

    Scalar radius() const
    { return radius_; }

    //- Tests
    bool isInside(const Point2D &testPoint) const;

    bool isOnEdge(const Point2D &point) const;

    bool isCovered(const Point2D &point) const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D &line) const;

    std::vector<Point2D> intersections(const LineSegment2D &line) const;

    Point2D nearestIntersect(const Point2D &point) const;

    LineSegment2D nearestEdge(const Point2D &point) const; // returns a unit tangent edge instead

    bool intersects(const Shape2D &shape) const;

    //- Transformations
    void scale(Scalar factor);

    void rotate(Scalar theta)
    {}

    //- Translations
    Circle &move(const Point2D& pos);

    Circle &operator+=(const Vector2D &translationVec);

    Circle &operator-=(const Vector2D &translationVec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

    //- Convenience
    Polygon polygonize(int nVerts) const;

    Polygon polygonize() const;

private:

    Point2D center_;
    Scalar radius_, area_, momentOfInertia_;

};

#endif
