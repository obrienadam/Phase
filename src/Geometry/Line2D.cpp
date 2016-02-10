#include "Line2D.h"

Line2D::Line2D(Scalar m, Scalar b)
    :
      m(m),
      b(b)
{

}

Line2D::Line2D(const Point2D &pt1, const Point2D &pt2)
{
    m = (pt2.y - pt1.y)/(pt2.x - pt1.x);
    b = pt1.y - m*pt1.x;
}

// External functions
Point2D lineIntersection(const Line2D &line1, const Line2D &line2)
{
    Scalar x = (line1.b - line2.b)/(line2.m - line1.m);
    return Point2D(x, line1(x));
}
