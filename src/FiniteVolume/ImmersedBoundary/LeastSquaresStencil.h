#ifndef LEAST_SQUARES_STENCIL_H
#define LEAST_SQUARES_STENCIL_H

#include "ImmersedBoundaryStencil.h"
#include "Matrix.h"

class LeastSquaresStencil: public ImmersedBoundaryStencil
{
public:
    LeastSquaresStencil(const Cell &cell, const Shape2D &shape, const CellGroup& cellGroup);

private:

    void constructMatrices();

    std::vector<Ref<const Cell>> cellPoints_;
    std::vector<Point2D> boundaryPoints_;
    std::vector<Point2D> boundaryNormals_;

    Matrix dMat_, nMat_;
};

#endif
