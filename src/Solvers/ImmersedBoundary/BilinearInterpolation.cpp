#include "BilinearInterpolation.h"

BilinearInterpolation::BilinearInterpolation(const FiniteVolumeGrid2D &grid)
{

}

BilinearInterpolation::BilinearInterpolation(const std::vector<Point2D>& pts)
{
    constructMatrix(pts);
}

Scalar BilinearInterpolation::operator()(const Scalar vals[], const Point2D& ip) const
{
    Matrix x(4, 1), phi(1, 4);

    x = {
      ip.x,
        ip.y,
        ip.x*ip.y,
        1.,
    };

    phi = {
      vals[0],
        vals[1],
        vals[2],
        vals[3],
    };

    return (phi*mat_*x)(0, 0);
}

std::vector<Scalar> BilinearInterpolation::operator()(const Point2D& ip) const
{
    Matrix x(4, 1), a(4, 1);

    x = {
      ip.x,
        ip.y,
        ip.x*ip.y,
        1.,
    };

    a = mat_*x;

    return a.containerCopy();
}

//- Private methods

void BilinearInterpolation::constructMatrix(const std::vector<Point2D> &pts)
{
    mat_.resize(pts.size(), 4);

    for(int i = 0; i < pts.size(); ++i)
    {
        mat_(i, 0) = pts[i].x;
        mat_(i, 1) = pts[i].y;
        mat_(i, 2) = pts[i].x*pts[i].y;
        mat_(i, 3) = 1.;
    }

    mat_.invert().transpose();
}
