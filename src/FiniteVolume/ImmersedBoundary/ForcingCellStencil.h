#ifndef FORCING_CELL_STENCIL_H
#define FORCING_CELL_STENCIL_H

#include "FiniteVolumeGrid2D.h"
#include "LinearInterpolation.h"
#include "ImmersedBoundaryStencil.h"

class ForcingCellStencil: public ImmersedBoundaryStencil
{
public:
    ForcingCellStencil(const Cell &cell, const Shape2D &shape, const CellGroup& cellGroup);

    const Cell& cell() const { return cell_; }

    const Point2D& xc() const { return xc_; }

    const std::vector<Ref<const Cell>>& iCells() const { return iCells_; }

    const std::vector<Scalar>& iCoeffs() const { return iCoeffs_; }

    const std::vector<Scalar>& nCoeffs() const { return nCoeffs_; }

    Scalar bCoeff() const { return bCoeff_; }

    Scalar bnCoeff() const { return bnCoeff_; }

private:

    std::vector<Ref<const Cell>> iCells_;
    std::vector<Scalar> iCoeffs_, nCoeffs_;
    Scalar bCoeff_, bnCoeff_;
    Point2D xc_;
    LinearInterpolation interpolator_;
};

#endif
