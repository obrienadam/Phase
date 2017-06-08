#ifndef IMMERSED_BOUNDARY_STENCIL_H
#define IMMERSED_BOUNDARY_STENCIL_H

#include "Cell.h"

class ImmersedBoundaryStencil
{
public:

    ImmersedBoundaryStencil(const Cell& cell) : cell_(cell) {}

    const Cell& cell() const { return cell_; }

protected:

    const Cell &cell_;
//    std::vector<Ref<const Cell>> stencilCells_;
//    std::vector<Point2D> boundaryPoints_;
//    std::vector<Scalar> cellCoeffs_;
//    std::vector<Scalar> boundaryCoeffs_;
};

#endif
