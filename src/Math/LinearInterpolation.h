#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H

#include "Interpolation.h"

class LinearInterpolation: public Interpolation
{
public:

    LinearInterpolation() {}
    LinearInterpolation(const std::vector<Point2D>& pts);

    Scalar operator()(const std::vector<Scalar> &vals, const Point2D &ip) const;

    Vector2D operator()(const std::vector<Vector2D> &vals, const Point2D &ip) const;

    std::vector<Scalar> operator()(const Point2D &ip) const;

    std::vector<Scalar> derivative(const Point2D &ip, const Vector2D& e, const std::vector<Point2D> &pts) const;

protected:

    virtual void constructMatrix(const std::vector<Point2D> &pts);
};

#endif
