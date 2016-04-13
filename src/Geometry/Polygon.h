#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include "Shape2D.h"
#include "Geometry.h"

class Polygon : public CGAL::Polygon_2<geometry::Kernel>
{
public:

    Polygon();
    Polygon(const std::vector<Point2D>& vertices);

    Point2D centroid() const { return centroid_; }
    Scalar area() const { return area_; }

    virtual bool isInside(const Point2D& testPoint) const;
    virtual bool isOnEdge(const Point2D& testPoint) const;

    virtual Point2D nearestIntersect(const Point2D& testPoint) const;

    bool isSimple() const;
    bool isConvex() const;

private:

    bool isSimple_, isConvex_;

    Point2D centroid_;
    Scalar area_;
};

#endif
