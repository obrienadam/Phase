#include "Shape2D.h"

std::vector<Point2D> Shape2D::intersections(const LineSegment2D &line) const
{
    std::vector<Point2D> points;

    for(const Point2D& point: intersections(Line2D(line.ptA(), line.norm())))
    {
        if(line.isBounded(point))
            points.push_back(point);
    }

    return points;
}