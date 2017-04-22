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

    template <class const_iterator>
    Polygon(const_iterator begin, const_iterator end)
    {
        for (const_iterator vtx = begin; vtx != end; ++vtx)
            boost::geometry::append(poly_, *vtx);

        init();
    }

    Polygon(const std::vector<Point2D> &vertices);

    Polygon(const boost::geometry::model::polygon<Point2D, false, true> &boostPgn);

    Type type() const { return POLYGON; }

    //- Polygon parameters
    const Point2D &centroid() const
    { return centroid_; }

    Scalar area() const
    { return area_; }

    const boost::geometry::model::polygon<Point2D, false, true> &boostPolygon() const
    { return poly_; }

    //- Tests
    virtual bool isInside(const Point2D &testPoint) const;

    virtual bool isOnEdge(const Point2D &testPoint) const;

    virtual bool isCovered(const Point2D &point) const;

    virtual bool isBoundedBy(const Point2D &point, Scalar toler) const;

    bool isValid() const;

    bool isEmpty() const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D &line) const;

    std::vector<Point2D> intersections(const LineSegment2D& line) const;

    Point2D nearestIntersect(const Point2D &point) const;

    LineSegment2D nearestEdge(const Point2D &point) const;

    bool intersects(const Shape2D &shape) const;

    //- Transformations
    void scale(Scalar factor);

    void rotate(Scalar theta);

    Polygon scale(Scalar factor) const;

    //- Translations
    Polygon &operator+=(const Vector2D &translationVec);

    Polygon &operator-=(const Vector2D &translationVec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

    //- Convenience
    Polygon polygonize() const
    { return *this; }

    //- Vertices
    const std::vector<Point2D> &vertices() const
    { return boost::geometry::exterior_ring(poly_); }

    std::vector<LineSegment2D> edges() const;

protected:

    void init();

    boost::geometry::model::polygon<Point2D, false, true> poly_;

    Point2D centroid_;
    Scalar area_;
};

//- External functions
Polygon intersectionPolygon(const Polygon &pgnA, const Polygon &pgnB);

Polygon difference(const Polygon &pgnA, const Polygon &pgnB);

Polygon clipPolygon(const Polygon &pgn, const Line2D &line);

#endif
