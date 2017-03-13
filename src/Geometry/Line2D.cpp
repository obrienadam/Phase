#include "Line2D.h"
#include "Matrix.h"
#include "Exception.h"

Line2D::Line2D()
{
    r0_ = n_ = d_ = Vector2D(0., 0.);
}

Line2D::Line2D(const Point2D &r0, const Vector2D &n)
        :
        r0_(r0),
        n_(n.unitVec())
{
    d_ = n_.tangentVec();
}

Point2D Line2D::operator()(Scalar t) const
{
    return r0_ + t * d_;
}

bool Line2D::isApproximatelyOnLine(const Point2D &pt) const
{
    Vector2D rVec = (pt - r0_).normalVec();
    return n_.isParallel(rVec);
}

bool Line2D::isAboveLine(const Point2D &pt) const
{
    return dot(pt - r0_, n_) > 0.;
}

bool Line2D::isBelowLine(const Point2D &pt) const
{
    return dot(pt - r0_, n_) < 0.;
}

//- Static functions

std::pair<Point2D, bool> Line2D::intersection(const Line2D &lineA, const Line2D &lineB)
{
    if (lineA.n().isParallel(lineB.n()))
        return std::make_pair(Point2D(), false); // lines are parallel

    Vector2D d1 = lineA.n().tangentVec(),
            d2 = lineB.n().tangentVec(),
            o1 = lineA.r0(),
            o2 = lineB.r0();

    Matrix A(2, 2), b(2, 1);

    A = {
            d1.x, -d2.x,
            d1.y, -d2.y,
    };

    b = {
            o2.x - o1.x,
            o2.y - o1.y,
    };

    A.solve(b);

    return std::make_pair(lineA(b(0, 0)), true);
}

//- External functions

std::ostream &operator<<(std::ostream &os, const Line2D &line)
{
    os << line.r0_ << " " << line.n_;
    return os;
}
