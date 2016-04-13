#include <math.h>

#include <CGAL/centroid.h>

#include "Polygon.h"

Polygon::Polygon()
{

}

Polygon::Polygon(const std::vector<Point2D> &vertices)
{
    typedef CGAL::Point_2<geometry::Kernel> Point;

    for(const Point2D& vtx: vertices)
        this->push_back(Point(vtx.x, vtx.y));

    isSimple_ = this->is_simple();
    isConvex_ = this->is_convex();

    area_ = CGAL::Polygon_2<geometry::Kernel>::area();

    Point centroid = CGAL::centroid(container().begin(), container().end(), CGAL::Dimension_tag<0>());
    centroid_ = Point2D(centroid.x(), centroid.y());
}

bool Polygon::isInside(const Point2D& testPoint) const
{
    return has_on_bounded_side(CGAL::Point_2<geometry::Kernel>(testPoint.x, testPoint.y));
}

bool Polygon::isOnEdge(const Point2D& testPoint) const
{
    return has_on_boundary(CGAL::Point_2<geometry::Kernel>(testPoint.x, testPoint.y));
}

Point2D Polygon::nearestIntersect(const Point2D& testPoint) const
{

}

bool Polygon::isSimple() const
{
    return isSimple_;
}

bool Polygon::isConvex() const
{
    return isConvex_;
}
