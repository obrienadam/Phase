#include <math.h>

#include <CGAL/centroid.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "Polygon.h"
#include "Exception.h"

Polygon::Polygon()
{

}

Polygon::Polygon(const std::vector<Point2D> &vertices)
{
    for(const Point2D& vtx: vertices)
        push_back(vtx.cgalPoint());

    init();
}

Polygon::Polygon(const CGAL::Polygon_2<Kernel> &other)
    :
      CGAL::Polygon_2<Kernel>::Polygon_2(other)
{
    init();
}

bool Polygon::isInside(const Point2D& testPoint) const
{
    return bounded_side(CGAL::Point_2<Kernel>(testPoint.cgalPoint())) == CGAL::ON_BOUNDED_SIDE;
}

bool Polygon::isOnEdge(const Point2D& testPoint) const
{
    return bounded_side(CGAL::Point_2<Kernel>(testPoint.cgalPoint())) == CGAL::ON_BOUNDARY;
}

Point2D Polygon::nearestIntersect(const Point2D& testPoint) const
{
    throw Exception("Polygon", "nearestIntersect", "not yet implemented.");
}

bool Polygon::isSimple() const
{
    return isSimple_;
}

bool Polygon::isConvex() const
{
    return isConvex_;
}

//- Protected methods

void Polygon::init()
{
    isSimple_ = is_simple();
    isConvex_ = is_convex();

    area_ = CGAL::Polygon_2<Kernel>::area();

    centroid_ = CGAL::centroid(container().begin(), container().end(), CGAL::Dimension_tag<0>());
}

//- External functions

bool doIntersect(const CGAL::Polygon_2<Kernel> &poly1, const CGAL::Polygon_2<Kernel> &poly2)
{
    return CGAL::do_intersect(poly1, poly2);
}

Polygon intersectionPolygon(const Polygon& poly1, const Polygon& poly2)
{
    std::vector< CGAL::Polygon_with_holes_2<Kernel> > result;
    CGAL::intersection(poly1, poly2, std::back_inserter(result));

    if(result.empty())
        throw Exception("", "intersectionPolygon", "polygons do not intersect.");
    else if(result.size() > 1)
        throw Exception("", "intersectionPolygon", "polygons intersect more than once.");

    return Polygon(result[0].outer_boundary());
}
