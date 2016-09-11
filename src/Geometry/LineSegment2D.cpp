#include "LineSegment2D.h"
#include "Line2D.h"

bool LineSegment2D::isBounded(const Point2D &pt) const
{
    return dot(pt - ptA_, ptB_ - ptA_) > 0 && (pt - ptA_).magSqr() <= magSqr();
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
