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

Point2D Line2D::intersection(const Line2D &other) const
{
    return (*this)(cross(other.r0_ - r0_, -other.d_)/cross(d_, -other.d_));
}

//- Static functions

//- External functions

std::ostream &operator<<(std::ostream &os, const Line2D &line)
{
    os << line.r0_ << " " << line.n_;
    return os;
}
