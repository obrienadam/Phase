#include "Quadrilateral.h"

Quadrilateral::Quadrilateral(const Point2D &v1, const Point2D &v2, const Point2D &v3, const Point2D &v4)
{
    vertices_[0] = v1;
    vertices_[1] = v2;
    vertices_[2] = v3;
    vertices_[3] = v4;
}

Scalar Quadrilateral::area() const
{
    Scalar area = 0.;

    for(int i = 1; i < 3; ++i)
        area += cross(vertices_[i] - vertices_[0], vertices_[i + 1] - vertices_[0]);

    return 0.5*area;
}
