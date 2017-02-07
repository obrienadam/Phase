#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Point2D.h"
#include "Matrix.h"

class Interpolation
{
public:

    virtual Scalar operator()(const std::vector<Scalar>& vals, const Point2D& ip) const = 0;
    virtual Vector2D operator()(const std::vector<Vector2D>& vals, const Point2D& ip) const = 0;

    virtual std::vector<Scalar> operator()(const Point2D& ip) const = 0;

    const Matrix& mat() const { return mat_; }

protected:

    virtual void constructMatrix(const std::vector<Point2D>& pts) = 0;

    Matrix mat_;

};

#include "BilinearInterpolation.h"
#include "QuadraticNormalInterpolation.h"

#endif
