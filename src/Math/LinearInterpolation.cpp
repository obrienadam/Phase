#include "LinearInterpolation.h"

LinearInterpolation::LinearInterpolation(const std::vector<Point2D> &pts)
{
    constructMatrix(pts);
}

Scalar LinearInterpolation::operator()(const std::vector<Scalar> &vals, const Point2D &ip) const
{
    Matrix b(1, 3);

    b = {
        vals[0],
        vals[1],
        vals[2]
    };

    Matrix x(3, 1);

    x = {
        ip.x,
        ip.y,
        1.
    };

    return (b*mat_*x)(0, 0);
}

Vector2D LinearInterpolation::operator()(const std::vector<Vector2D> &vals, const Point2D &ip) const
{
    Matrix bx(1, 3), by(1, 3);

    bx = {
        vals[0].x,
        vals[1].x,
        vals[2].x
    };

    by = {
        vals[0].y,
        vals[1].y,
        vals[2].y
    };

    Matrix x(3, 1);

    x = {
        ip.x,
        ip.y,
        1.
    };

    return Vector2D((mat_*bx*x)(0, 0), (mat_*by*x)(0, 0));
}

std::vector<Scalar> LinearInterpolation::operator()(const Point2D &ip) const
{
    Matrix x(3, 1), a(3, 1);

    x = {
            ip.x,
            ip.y,
            1.,
    };

    a = mat_ * x;

    return a.containerCopy();
}

std::vector<Scalar> LinearInterpolation::derivative(const Point2D& ip, const Vector2D &e, const std::vector<Point2D>& pts) const
{
    Matrix mat(3, 3), x(3, 1), a(3, 1);

    mat(0, 0) = e.x;
    mat(0, 1) = e.y;
    mat(0, 2) = 0.;

    for(int i = 1; i < 3; ++i)
    {
        mat(i, 0) = pts[i].x;
        mat(i, 1) = pts[i].y;
        mat(i, 2) = 1.;
    }

    mat.invert().transpose();

    x = {
        ip.x,
        ip.y,
        1.
    };

    a = mat*x;

    return a.containerCopy();
}

void LinearInterpolation::constructMatrix(const std::vector<Point2D> &pts)
{
    mat_.resize(3, 3);

    for (int i = 0; i < 3; ++i)
    {
        mat_(i, 0) = pts[i].x;
        mat_(i, 1) = pts[i].y;
        mat_(i, 2) = 1.;
    }

    mat_.invert().transpose();
}
