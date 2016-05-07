#include "Line2D.h"
#include "Matrix.h"
#include "Exception.h"

Line2D::Line2D(const Point2D &r0, const Point2D &d)
    :
      r0_(r0),
      d_(d.unitVec())
{
    n_ = d_.normalVec();
}

Point2D Line2D::operator()(Scalar t) const
{
    return r0_ + t*d_;
}

bool Line2D::isApproximatelyOnLine(const Point2D& pt) const
{
    Scalar t = (pt.x - r0_.x)/d_.x;
    return fabs(r0_.y + t*d_.y - pt.y) < 1e-10;
}

bool Line2D::isAboveLine(const Point2D& pt) const
{
    return dot(pt - r0_, n_) > 0.;
}

bool Line2D::isBelowLine(const Point2D& pt) const
{
    return dot(pt - r0_, n_) < 0.;
}

Scalar Line2D::t(const Point2D &pt) const
{
    if(isApproximatelyOnLine(pt))
        return 0.5*((pt.x - r0_.x)/d_.x + (pt.y - r0_.y)/d_.y);
    else
        throw Exception("Line2D", "t", "specified point is not on the line.");
}

//- Static functions

Point2D Line2D::intersection(const Line2D &lineA, const Line2D &lineB)
{
    Matrix A(2, 2), b(2, 1);

    A(0, 0) = lineA.d_.x;
    A(0, 1) = -lineB.d_.x;
    A(1, 0) = lineA.d_.y;
    A(1, 1) = -lineB.d_.y;

    b(0, 0) = lineB.r0_.x - lineA.r0_.x;
    b(1, 0) = lineB.r0_.y - lineA.r0_.y;

    A.solve(b);

    return lineA(b(0, 0));
}

//- External functions

std::ostream& operator<<(std::ostream& os, const Line2D& line)
{
    os << line.r0_ << " " << line.d_;
    return os;
}
