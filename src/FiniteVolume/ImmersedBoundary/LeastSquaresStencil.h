#ifndef LEAST_SQUARES_STENCIL_H
#define LEAST_SQUARES_STENCIL_H

#include "ImmersedBoundaryStencil.h"
#include "Matrix.h"

class LeastSquaresStencil: public ImmersedBoundaryStencil
{
public:
    LeastSquaresStencil(const Cell &cell, const Shape2D &shape);

    const std::vector<Ref<const Cell>>& iCells() const { return cellPoints_; }

    const std::vector<std::pair<Point2D, Vector2D>>& boundaryPoints() const { return boundaryPoints_; }

    const std::vector<Scalar>& dirichletCellCoeffs() const { return dirichletCellCoeffs_; }

    const std::vector<Scalar>& dirichletBoundaryCoeffs() const { return dirichletBoundaryCoeffs_; }

    const std::vector<Scalar>& neumannCellCoeffs() const { return neumannCellCoeffs_; }

    const std::vector<Scalar>& neumannBoundaryCoeffs() const { return neumannBoundaryCoeffs_; }

private:

    std::vector<Ref<const Cell>> cellPoints_;
    std::vector<std::pair<Point2D, Vector2D>> boundaryPoints_; //- contains a point and a unit normal

    std::vector<Scalar> dirichletCellCoeffs_, neumannCellCoeffs_;
    std::vector<Scalar> dirichletBoundaryCoeffs_, neumannBoundaryCoeffs_;

    Point2D xc_;
    Matrix dMat_, nMat_;
};

#endif
