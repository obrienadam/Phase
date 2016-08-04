#ifndef CIRCLE_H
#define CIRCLE_H

#include "Types.h"
#include "Shape2D.h"
#include "Polygon.h"

class Circle : public Shape2D
{
public:

    Circle(Point2D center, Scalar radius);

    void init(const Point2D &center, Scalar radius);

    const Point2D& centroid() const { return center_; }
    Scalar area() const;
    Scalar radius() const { return radius_; }

    void scale(Scalar factor);
    void rotate(Scalar theta) {}

    bool isInside(const Point2D &testPoint) const;
    bool isOnEdge(const Point2D &testPoint) const;
    Point2D nearestIntersect(const Point2D &testPoint) const;
    std::pair<Point2D, bool> firstIntersect(Point2D ptA, Point2D ptB) const;

    boost::geometry::model::box<Point2D> boundingBox() const;

    void operator+=(const Vector2D& translationVec);
    void operator-=(const Vector2D& translationVec);

private:

    Point2D center_;
    Scalar radius_, area_;

};

#endif
