#ifndef QUADRATIC_NORMAL_INTERPOLATION_H
#define QUADRATIC_NORMAL_INTERPOLATION_H

#include "Interpolation.h"

class QuadraticNormalInterpolation : public Interpolation
{
public:
    QuadraticNormalInterpolation(){}
    QuadraticNormalInterpolation(const std::vector<Point2D> &pts, const Vector2D& normal);

    Scalar operator()(const std::vector<Scalar>& vals, const Point2D& ip) const;
    Vector2D operator()(const std::vector<Vector2D>& vals, const Point2D& ip) const;

    std::vector<Scalar> operator()(const Point2D& ip) const;

private:

    void constructMatrix(const std::vector<Point2D>& pts);

    Vector2D normal_;
};

#endif
