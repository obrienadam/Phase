#include "Line2D.h"
#include "Matrix.h"
#include "Exception.h"

Line2D::Line2D(const Point2D &r0, const Point2D &n)
    :
      r0_(r0),
      n_(n.unitVec())
{
    d_ = n_.tangentVec();
}

Point2D Line2D::operator()(Scalar t) const
{
    return r0_ + t*d_;
}

bool Line2D::isApproximatelyOnLine(const Point2D& pt) const
{
    Vector2D rVec = (pt - r0_).normalVec();
    Scalar proj = dot(rVec, n_);

    return fabs(proj*proj - rVec.magSqr()) < std::numeric_limits<Scalar>::epsilon() ? true : false;
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
    Scalar proj = dot(lineA.n(), lineB.n());

    if(fabs(proj*proj - 1) < std::numeric_limits<Scalar>::epsilon())
        return std::make_pair(Point2D(), false); // lines are parallel

    Vector2D d1 = lineA.n().tangentVec(),
            d2 = lineB.n().tangentVec(),
            o1 = lineA.r0(),
            o2 = lineB.r0();

    Matrix A(2, 2), b(2, 1);

    A(0, 0) = d1.x;
    A(0, 1) = -d2.x;
    A(1, 0) = d1.y;
    A(1, 1) = -d2.y;

    b(0, 0) = o2.x - o1.x;
    b(1, 0) = o2.y - o1.y;

    A.solve(b);

    return std::make_pair(lineA(b(0, 0)), true);
}

//- External functions

std::ostream& operator<<(std::ostream& os, const Line2D& line)
{
    os << line.r0_ << " " << line.n_;
    return os;
}
