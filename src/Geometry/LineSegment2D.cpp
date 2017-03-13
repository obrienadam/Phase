#include "LineSegment2D.h"

bool LineSegment2D::isBounded(const Point2D &pt) const
{
    return dot(pt - ptA_, ptB_ - ptA_) > 0 - Point2D::epsilon() &&
           (pt - ptA_).magSqr() <= lengthSqr() + Point2D::epsilon();
}
