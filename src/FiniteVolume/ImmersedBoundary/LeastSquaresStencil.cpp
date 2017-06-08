#include "LeastSquaresStencil.h"
#include "CellGroup.h"
#include "Exception.h"

LeastSquaresStencil::LeastSquaresStencil(const Cell &cell, const Shape2D &shape)
    :
      ImmersedBoundaryStencil(cell)
{
    Vector2D xc, n;

    xc = shape.nearestIntersect(cell.centroid());
    n = (cell.centroid() - xc).unitVec();

    xc_ = xc;

    boundaryPoints_.push_back(std::make_pair(xc, n));

    auto hasNeighbourInSolid = [&shape](const Cell& cell) {
        for(const InteriorLink& nb: cell.neighbours())
            if(shape.isInside(nb.cell().centroid()))
                return true;

        return false;
    };

    for(const InteriorLink& nb: cell.neighbours())
        if(!shape.isInside(nb.cell().centroid())) // make sure it is not a solid cell
        {
            cellPoints_.push_back(nb.cell());

            if(hasNeighbourInSolid(nb.cell()))
            {
                xc = shape.nearestIntersect(nb.cell().centroid());
                n = (nb.cell().centroid() - xc).unitVec();

                boundaryPoints_.push_back(std::make_pair(xc, n));
            }
        }

    for(const DiagonalCellLink& dg: cell.diagonals())
        if(!shape.isInside(dg.cell().centroid()))
        {
            cellPoints_.push_back(dg.cell());

            if(hasNeighbourInSolid(dg.cell()))
            {
                xc = shape.nearestIntersect(dg.cell().centroid());
                n = (dg.cell().centroid() - xc).unitVec();

                boundaryPoints_.push_back(std::make_pair(xc, n));
            }
        }

    if(cellPoints_.size() + boundaryPoints_.size() < 6)
        throw Exception("LeastSquaresStencil", "LeastSquaresStencil", "stencil is too small.");

    dMat_ = Matrix(cellPoints_.size() + boundaryPoints_.size(), 6);
    nMat_ = Matrix(cellPoints_.size() + boundaryPoints_.size(), 6);

    int i = 0;
    for(const Cell& cell: cellPoints_)
    {
        Point2D pt = cell.centroid();

        dMat_(i, 0) = pt.x*pt.x;
        dMat_(i, 1) = pt.y*pt.y;
        dMat_(i, 2) = pt.x*pt.y;
        dMat_(i, 3) = pt.x;
        dMat_(i, 4) = pt.y;
        dMat_(i, 5) = 1.;

        nMat_(i, 0) = pt.x*pt.x;
        nMat_(i, 1) = pt.y*pt.y;
        nMat_(i, 2) = pt.x*pt.y;
        nMat_(i, 3) = pt.x;
        nMat_(i, 4) = pt.y;
        nMat_(i, 5) = 1.;

        i++;
    }

    for(const auto& bp: boundaryPoints_)
    {
        xc = bp.first;
        n = bp.second;

        dMat_(i, 0) = xc.x*xc.x;
        dMat_(i, 1) = xc.y*xc.y;
        dMat_(i, 2) = xc.x*xc.y;
        dMat_(i, 3) = xc.x;
        dMat_(i, 4) = xc.y;
        dMat_(i, 5) = 1.;

        nMat_(i, 0) = 2*xc.x*n.x;
        nMat_(i, 1) = 2*xc.y*n.y;
        nMat_(i, 2) = xc.y*n.x + xc.x*n.y;
        nMat_(i, 3) = n.x;
        nMat_(i, 4) = n.y;
        nMat_(i, 5) = 0.;

        i++;
    }

    dMat_ = dMat_*transpose(inverse(transpose(dMat_)*dMat_));
    nMat_ = nMat_*transpose(inverse(transpose(nMat_)*nMat_));

    Matrix x(6, 1);

    Point2D pt = cell.centroid();

    x = {
        pt.x*pt.x,
        pt.y*pt.y,
        pt.x*pt.y,
        pt.x,
        pt.y,
        1.
    };

    std::vector<Scalar> coeffs = (dMat_*x).containerCopy();

    dirichletCellCoeffs_.insert(dirichletCellCoeffs_.begin(), coeffs.begin(), coeffs.begin() + cellPoints_.size());
    dirichletBoundaryCoeffs_.insert(dirichletBoundaryCoeffs_.begin(), coeffs.begin() + cellPoints_.size(), coeffs.end());

    coeffs = (nMat_*x).containerCopy();

    neumannCellCoeffs_.insert(neumannCellCoeffs_.begin(), coeffs.begin(), coeffs.begin() + cellPoints_.size());
    neumannBoundaryCoeffs_.insert(neumannBoundaryCoeffs_.begin(), coeffs.begin() + cellPoints_.size(), coeffs.end());
}
