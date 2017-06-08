#ifndef FORCING_CELL_STENCIL_H
#define FORCING_CELL_STENCIL_H

#include "FiniteVolumeGrid2D.h"
#include "ImmersedBoundaryStencil.h"
#include "Matrix.h"

class ForcingCellStencil: public ImmersedBoundaryStencil
{
public:
    ForcingCellStencil(const Cell &cell, const Shape2D &shape, const CellGroup& cellGroup);

    const Point2D& xc() const { return xc_; }

    const std::vector<Ref<const Cell>>& nbCells() const { return nbCells_; }

    const std::vector<Scalar>& dirichletCellCoeffs() const { return dirichletCellCoeffs_; }

    const std::vector<Scalar>& neumannCellCoeffs() const { return neumannCellCoeffs_; }

    Scalar dirichletBoundaryCoeff() const { return dirichletBoundaryCoeff_; }

    Scalar neumannBoundaryCoeff() const { return neumannBoundaryCoeff_; }

private:

    std::vector<Ref<const Cell>> nbCells_;
    std::vector<Scalar> dirichletCellCoeffs_, neumannCellCoeffs_;
    Scalar dirichletBoundaryCoeff_, neumannBoundaryCoeff_;
    Point2D xc_;
};

#endif
