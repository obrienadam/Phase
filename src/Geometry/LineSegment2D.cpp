#include "LineSegment2D.h"
#include "Line2D.h"

bool LineSegment2D::isBounded(const Point2D &pt) const
{
    return dot(pt - ptA_, ptB_ - ptA_) > 0 - Point2D::epsilon() && (pt - ptA_).magSqr() <= magSqr() + Point2D::epsilon();
}

//- External functions

std::vector<Point2D> intersection(const LineSegment2D& lineSeg, const Shape2D& shape)
{
    Line2D line(lineSeg.ptA(), lineSeg.norm());

    std::vector<Point2D> xc = shape.intersections(line);
    std::vector<Point2D> result;

    for(const Point2D& pt: xc)
    {
        if(lineSeg.isBounded(pt))
            result.push_back(pt);
    }

    return result;
}

Point2D nearest(const LineSegment2D& lineSeg, const Shape2D& shape)
{
    auto xc = intersection(lineSeg, shape);

    if(xc.size() != 0)
        return xc[0];

    Point2D pt1 = shape.nearestIntersect(lineSeg.ptA()), pt2 = shape.nearestIntersect(lineSeg.ptB());

    return (lineSeg.ptA() - pt1).magSqr() < (lineSeg.ptB() - pt2).magSqr() ? lineSeg.ptA() : lineSeg.ptB();
}
