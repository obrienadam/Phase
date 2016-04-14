#include "LineSegment2D.h"

LineSegment2D::LineSegment2D(const Point2D &pt1, const Point2D &pt2)
    :
      CGAL::Segment_2<geometry::Kernel>(
          CGAL::Point_2<geometry::Kernel>(pt1.x, pt1.y),
          CGAL::Point_2<geometry::Kernel>(pt2.x, pt2.y)
          )
{

}

//- External functions

std::pair<Point2D, bool> intersection(LineSegment2D &seg1, const LineSegment2D &seg2)
{
    CGAL::Object result = CGAL::intersection(seg1, seg2);

    if(const CGAL::Point_2<geometry::Kernel> *pt = CGAL::object_cast<CGAL::Point_2<geometry::Kernel> >(&result))
        return std::pair<Point2D, bool>(Point2D(pt->x(), pt->y()), true);
    else
        return std::pair<Point2D, bool>(Point2D(), false);
}
