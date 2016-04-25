#include "BilinearInterpolation.h"

BilinearInterpolation::BilinearInterpolation(const FiniteVolumeGrid2D &grid)
{
    // Finite the four nearest neighbour cells
    Point2D centroids[4];

    // Construct
    constructMatrix(centroids);
}

BilinearInterpolation::BilinearInterpolation(const Point2D pts[])
{
    constructMatrix(pts);
}

Scalar BilinearInterpolation::operator()(const Scalar vals[], const Point2D& ip) const
{
    Matrix x(4, 1), phi(1, 4);

    x(0, 0) = ip.x;
    x(1, 0) = ip.y;
    x(2, 0) = ip.x*ip.y;
    x(3, 0) = 1.;

    phi.init(vals, vals + 4);

    return (phi*mat_*x)(0, 0);
}

std::vector<Scalar> BilinearInterpolation::operator()(const Point2D& ip) const
{
    Matrix x(4, 1), a(4, 1);

    x(0, 0) = ip.x;
    x(1, 0) = ip.y;
    x(2, 0) = ip.x*ip.y;
    x(3, 0) = 1.;

    a = mat_*x;

    return a.containerCopy();
}

//- Private methods

void BilinearInterpolation::constructMatrix(const Point2D pts[])
{
    mat_.resize(4, 4);

    for(int i = 0; i < 4; ++i)
    {
        mat_(i, 0) = pts[i].x;
        mat_(i, 1) = pts[i].y;
        mat_(i, 2) = pts[i].x*pts[i].y;
        mat_(i, 3) = 1.;
    }
    mat_.invert().transpose();
}
