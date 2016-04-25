#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <CGAL/Polygon_2.h>

#include "Shape2D.h"
#include "Geometry.h"

class Polygon : public CGAL::Polygon_2<Kernel>
{
public:

    Polygon();
    Polygon(const std::vector<Point2D>& vertices);
    Polygon(const CGAL::Polygon_2<Kernel>& other);

    template<class T>
    Polygon(T begin, T end);

    const Point2D& centroid() const { return centroid_; }
    Scalar area() const { return area_; }

    virtual bool isInside(const Point2D& testPoint) const;
    virtual bool isOnEdge(const Point2D& testPoint) const;

    virtual Point2D nearestIntersect(const Point2D& testPoint) const;

    bool isSimple() const;
    bool isConvex() const;

protected:

    void init();

    bool isSimple_, isConvex_;

    Point2D centroid_;
    Scalar area_;
};

bool doIntersect(const CGAL::Polygon_2<Kernel> &poly1, const CGAL::Polygon_2<Kernel> &poly2);
Polygon intersectionPolygon(const Polygon& poly1, const Polygon& poly2);

#endif
