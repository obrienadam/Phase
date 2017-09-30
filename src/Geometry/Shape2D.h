#ifndef SHAPE_2D_H
#define SHAPE_2D_H

#include "Point2D.h"
#include "Line2D.h"
#include "LineSegment2D.h"
#include "Ray2D.h"

class Polygon;

class Shape2D
{
public:

    enum Type{POLYGON, CIRCLE, BOX};

    virtual ~Shape2D()
    {}

    //- Type
    virtual Type type() const = 0;

    //- Shape parameters
    virtual const Point2D &centroid() const = 0;

    virtual Scalar area() const = 0;

    virtual Scalar momentOfInertia() const = 0;

    virtual Scalar perimeter() const = 0;

    //- Tests
    virtual bool isInside(const Point2D &point) const = 0;

    virtual bool isOutside(const Point2D &point) const
    { return !isInside(point); }

    virtual bool isOnEdge(const Point2D &point) const = 0;

    virtual bool isCovered(const Point2D &point) const = 0;

    template <class const_iterator>
    Point2D closest(const_iterator begin, const_iterator end) const
    {
        Point2D pt;
        Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

        for(const_iterator itr = begin; itr != end; ++itr)
        {
            Scalar distSqr = (nearestIntersect(*itr) - *itr).magSqr();

            if(distSqr < minDistSqr)
                pt = *itr;
        }

        return pt;
    }

    //- Intersections
    virtual std::vector<Point2D> intersections(const Line2D &line) const = 0;

    virtual std::vector<Point2D> intersections(const LineSegment2D &line) const = 0;

    virtual std::vector<Point2D> intersections(const Ray2D& ray) const = 0;

    virtual Point2D nearestIntersect(const Point2D &point) const = 0;

    virtual Point2D nearestPoint(const Shape2D& shape) const;

    virtual LineSegment2D nearestEdge(const Point2D &point) const = 0;

    virtual bool intersects(const Shape2D &shape) const = 0;

    //- Transformations
    virtual void scale(Scalar factor) = 0;

    virtual void rotate(Scalar theta) = 0;

    //- Translations
    virtual Shape2D &move(const Point2D& pos) = 0;

    virtual Shape2D &operator+=(const Vector2D &vec) = 0;

    virtual Shape2D &operator-=(const Vector2D &vec) = 0;

    //- Bounding box
    virtual boost::geometry::model::box<Point2D> boundingBox() const = 0;

    //- Convenience
    virtual Polygon polygonize() const = 0;

private:
};

#endif
