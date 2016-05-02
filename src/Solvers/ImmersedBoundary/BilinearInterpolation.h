#ifndef BILINEAR_INTERPOLATION_H
#define BILINEAR_INTERPOLATION_H

#include "Matrix.h"
#include "FiniteVolumeGrid2D.h"

class BilinearInterpolation
{
public:

    BilinearInterpolation(const FiniteVolumeGrid2D& grid);
    BilinearInterpolation(const std::vector<Point2D> &pts);

    Scalar operator()(const Scalar vals[], const Point2D& ip) const;
    std::vector<Scalar> operator()(const Point2D& ip) const;

    const Matrix& mat() const { return mat_; }

private:

    void constructMatrix(const std::vector<Point2D>& pts);

    Matrix mat_;

};

#endif
