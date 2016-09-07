#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <boost/geometry/geometries/polygon.hpp>

#include "Shape2D.h"
#include "Line2D.h"

class Polygon : public Shape2D
{
public:

    Polygon();
    Polygon(const std::vector<Point2D>& vertices);
    Polygon(const boost::geometry::model::polygon<Point2D, false, true>& boostPgn);

    //- Polygon parameters
    const Point2D& centroid() const { return centroid_; }
    Scalar area() const { return area_; }
    const boost::geometry::model::polygon<Point2D, false, true>& boostPolygon() const { return poly_; }

    //- Tests
    virtual bool isInside(const Point2D& testPoint) const;
    virtual bool isOnEdge(const Point2D& testPoint) const;
    virtual bool isCovered(const Point2D &point) const;

    bool isValid() const;
    bool isEmpty() const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D& line) const;
    Point2D nearestIntersect(const Point2D &point) const;
    std::pair<Point2D, Point2D> nearestEdge(const Point2D& point) const;

    //- Transformations
    void scale(Scalar factor);
    void rotate(Scalar theta);
    Polygon scale(Scalar factor) const;

    //- Translations
    virtual void operator+=(const Vector2D& translationVec);
    virtual void operator-=(const Vector2D& translationVec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

    //- Convenience
    Polygon polygonize() const { return *this; }

    //- Vertices
    const std::vector<Point2D>& vertices() const { return boost::geometry::exterior_ring(poly_); }

protected:

    void init();

    boost::geometry::model::polygon<Point2D, false, true> poly_;

    Point2D centroid_;
    Scalar area_;
};

//- External functions
Polygon intersectionPolygon(const Polygon &pgnA, const Polygon &pgnB);
Polygon difference(const Polygon &pgnA, const Polygon &pgnB);
Polygon clipPolygon(const Polygon& pgn, const Line2D& line);

#endif
