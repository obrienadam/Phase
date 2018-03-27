#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <boost/geometry/geometries/polygon.hpp>

#include "System/Exception.h"

#include "Shape2D.h"
#include "Line2D.h"

class Polygon : public Shape2D
{
public:

    template <class const_iterator>
    static Polygon convexHull(const_iterator begin, const_iterator end)
    {
        boost::geometry::model::ring<Point2D, false, true> poly;
        boost::geometry::model::multi_point<Point2D> pts(begin, end);
        boost::geometry::convex_hull(pts, poly);
        return Polygon(poly);
    }

    Polygon();

    template<class const_iterator>
    Polygon(const_iterator begin, const_iterator end)
    {
        for (const_iterator vtx = begin; vtx != end; ++vtx)
            boost::geometry::append(poly_, *vtx);

        init();
    }

    Polygon(const std::initializer_list<Point2D> &vertices);

    Polygon(const boost::geometry::model::ring<Point2D, false, true> &boostRing);

    Type type() const
    { return POLYGON; }

    //- Polygon parameters
    const Point2D &centroid() const
    { return centroid_; }

    Scalar area() const
    { return area_; }

    Scalar momentOfInertia() const
    { throw Exception("Polygon", "momentOfInertia", "not implemented."); }

    Scalar perimeter() const;

    const boost::geometry::model::ring<Point2D, false, true> &boostRing() const
    { return poly_; }

    //- Tests
    virtual bool isInside(const Point2D &testPoint) const;

    virtual bool isOnEdge(const Point2D &testPoint) const;

    virtual bool isCovered(const Point2D &point) const;

    bool isValid() const;

    bool isEmpty() const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D &line) const;

    std::vector<Point2D> intersections(const LineSegment2D &line) const;

    std::vector<Point2D> intersections(const Ray2D& ray) const;

    Point2D nearestIntersect(const Point2D &point) const;

    LineSegment2D nearestEdge(const Point2D &point) const;

    bool intersects(const Shape2D &shape) const;

    //- Transformations
    void scale(Scalar factor);

    void rotate(Scalar theta);

    Polygon scale(Scalar factor) const;

    //- Translations
    Polygon &move(const Point2D& pos);

    Polygon &operator+=(const Vector2D &translationVec);

    Polygon &operator-=(const Vector2D &translationVec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

    //- Convenience
    Polygon polygonize() const
    { return *this; }

    //- Vertices
    const std::vector<Point2D>& vertices() const
    { return poly_; }

    std::vector<LineSegment2D> edges() const;

protected:

    void init();

    boost::geometry::model::ring<Point2D, false, true> poly_;

    Point2D centroid_;
    Scalar area_;

    bool valid_, simple_;
};

//- External functions
std::vector<Polygon> intersection(const Polygon &pgnA, const Polygon &pgnB);

std::vector<Polygon> difference(const Polygon &pgnA, const Polygon &pgnB);

#endif
