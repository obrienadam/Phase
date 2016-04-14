#ifndef LINE_SEGMENT_2D_H
#define LINE_SEGMENT_2D_H

#include "Geometry.h"
#include "Point2D.h"

class LineSegment2D: public CGAL::Segment_2<geometry::Kernel>
{
public:
    LineSegment2D(const Point2D& pt1, const Point2D& pt2);
};

std::pair<Point2D, bool> intersection(LineSegment2D &seg1, const LineSegment2D &seg2);

#endif
