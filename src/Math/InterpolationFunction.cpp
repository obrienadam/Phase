#include "InterpolationFunction.h"

InterpolationFunction::InterpolationFunction(const std::vector<Point2D> &pts)
{
    mat_ = Matrix(pts.size(), pts.size());

    int i = 0;
    for(const Point2D& pt: pts)
    {
        int j = 0;
        mat_(i, j++) = pt.x*pt.y*pt.y;
        mat_(i, j++) = pt.y*pt.y;
        mat_(i, j++) = pt.x*pt.y;
        mat_(i, j++) = pt.x;
        mat_(i, j++) = pt.y;
        mat_(i++, j++) = 1.;
    }

    mat_.invert();
}

InterpolationFunction::InterpolationFunction(const std::vector<Point2D> &pts, const std::vector<Scalar> &vals)
    :
      InterpolationFunction(pts)
{
    b_ = Matrix(vals.size(), 1);
    c_ = mat_*b_;
}

Scalar InterpolationFunction::operator()(const Point2D &pt) const
{
    Matrix x(6, 1);

    x(0, 0) = pt.x*pt.y*pt.y;
    x(1, 0) = pt.y*pt.y;
    x(2, 0) = pt.x*pt.y;
    x(3, 0) = pt.x;
    x(4, 0) = pt.y;
    x(5, 0) = 1.;

    return (transpose(c_)*x)(0, 0);
}
