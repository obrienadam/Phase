#ifndef LINE_SEGMENT_2D_H
#define LINE_SEGMENT_2D_H

#include "Point2D.h"
#include "Shape2D.h"

class LineSegment2D
{
public:

    LineSegment2D(const Point2D& ptA, const Point2D& ptB) : ptA_(ptA), ptB_(ptB) {}

    const Point2D& ptA() const { return ptA_; }
    const Point2D& ptB() const { return ptB_; }

    Scalar magSqr() const { return (ptB_ - ptA_).magSqr(); }
    Scalar mag() const { return (ptB_ - ptA_).mag(); }

    Point2D center() const { return (ptA_ + ptB_)/2.; }
    Vector2D norm() const { return (ptB_ - ptA_).normalVec(); }

    bool isBounded(const Point2D& pt) const;

private:

    Point2D ptA_, ptB_;

};

std::vector<Point2D> intersection(const LineSegment2D& lineSeg, const Shape2D& shape);
Point2D nearest(const LineSegment2D& lineSeg, const Shape2D& shape);

#endif
