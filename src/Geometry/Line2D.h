#ifndef LINE_2D_H
#define LINE_2D_H

#include "Types.h"
#include "Point2D.h"

class Line2D
{
public:

    Line2D(Scalar m = 0., Scalar b = 0.);
    Line2D(const Point2D& pt1, const Point2D& pt2);

    Scalar operator()(Scalar x) const { return m*x + b; }

    Scalar m, b;
};

Point2D lineIntersection(const Line2D& line1, const Line2D& line2);

#endif
