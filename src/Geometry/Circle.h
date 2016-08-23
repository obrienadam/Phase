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

    void init(const Point2D &center, Scalar radius);

    //- Circle parameters
    const Point2D& centroid() const { return center_; }
    Scalar area() const { return area_; }
    Scalar radius() const { return radius_; }

    //- Tests
    bool isInside(const Point2D &testPoint) const;
    bool isOnEdge(const Point2D &testPoint) const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D &line) const;
    Point2D nearestIntersect(const Point2D &point) const;
    std::pair<Point2D, Point2D> nearestEdge(const Point2D &point) const; // returns a unit tangent edge instead

    //- Transformations
    void scale(Scalar factor);
    void rotate(Scalar theta) { }

    //- Translations
    void operator+=(const Vector2D& translationVec);
    void operator-=(const Vector2D& translationVec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

private:

    Point2D center_;
    Scalar radius_, area_;

};

#endif
