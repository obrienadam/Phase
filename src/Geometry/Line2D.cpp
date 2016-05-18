#include "Line2D.h"
#include "Matrix.h"
#include "Exception.h"

Line2D::Line2D(const Point2D &r0, const Point2D &n)
    :
      r0_(r0),
      n_(n.unitVec())
{
    d_ = n_.normalVec();
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

//- Static functions

std::pair<Point2D, bool> Line2D::intersection(const Line2D &lineA, const Line2D &lineB)
{
    Vector2D da = lineA.n().abs(),
            db = lineB.n().abs();

    const double toler = 1e-12;
    if(fabs(da.x - db.x) < toler
            && fabs(da.y - db.y) < toler)
        return std::make_pair(Point2D(INFINITY, INFINITY), false);

    Matrix A(2, 2), b(2, 1);
    da = lineA.d_;
    db = lineB.d_;

    A(0, 0) = da.x;
    A(0, 1) = -db.x;
    A(1, 0) = da.y;
    A(1, 1) = -db.y;

    b(0, 0) = lineB.r0_.x - lineA.r0_.x;
    b(1, 0) = lineB.r0_.y - lineA.r0_.y;

    A.solve(b);

    return std::make_pair(lineA(b(0, 0)), true);
}

//- External functions

std::ostream& operator<<(std::ostream& os, const Line2D& line)
{
    os << line.r0_ << " " << line.n_;
    return os;
}
