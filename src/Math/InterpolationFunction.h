#ifndef INTERPOLATION_FUNCTION_H
#define INTERPOLATION_FUNCTION_H

#include "Matrix.h"
#include "Point2D.h"

class InterpolationFunction
{
public:
    InterpolationFunction(const std::vector<Point2D>& pts);
    InterpolationFunction(const std::vector<Point2D>& pts, const std::vector<Scalar>& vals);

    Scalar operator()(const Point2D& pt) const;

private:

    Matrix mat_, c_, b_;
};

#endif
