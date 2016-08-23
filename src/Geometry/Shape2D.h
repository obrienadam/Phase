#ifndef SHAPE_2D_H
#define SHAPE_2D_H

#include "Point2D.h"
#include "Line2D.h"

class Shape2D
{
public:

    virtual ~Shape2D() {}

    //- Shape parameters
    virtual const Point2D& centroid() const = 0;
    virtual Scalar area() const = 0;

    //- Tests
    virtual bool isInside(const Point2D& point) const = 0;
    virtual bool isOutside(const Point2D& point) const { return !isInside(point); }
    virtual bool isOnEdge(const Point2D& point) const = 0;

    //- Intersections
    virtual std::vector<Point2D> intersections(const Line2D& line) const = 0;
    virtual Point2D nearestIntersect(const Point2D& point) const = 0;
    virtual std::pair<Point2D, Point2D> nearestEdge(const Point2D& point) const = 0;

    //- Transformations
    virtual void scale(Scalar factor) = 0;
    virtual void rotate(Scalar theta) = 0;

    //- Translations
    virtual void operator+=(const Vector2D& vec) = 0;
    virtual void operator-=(const Vector2D& vec) = 0;

    //- Bounding box
    virtual boost::geometry::model::box<Point2D> boundingBox() const = 0;

private:
};

#endif
