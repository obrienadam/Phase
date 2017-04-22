#include "LineSegment2D.h"

bool LineSegment2D::isBounded(const Point2D &pt) const
{
    return dot(pt - ptA_, ptB_ - ptA_) >= 0. && dot(pt - ptA_, ptB_ - ptA_) <= (ptB_ - ptA_).magSqr();
}
