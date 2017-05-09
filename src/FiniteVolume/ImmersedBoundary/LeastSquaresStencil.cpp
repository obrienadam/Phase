#include "LeastSquaresStencil.h"
#include "CellGroup.h"

LeastSquaresStencil::LeastSquaresStencil(const Cell &cell, const Shape2D &shape, const CellGroup &cellGroup)
    :
      ImmersedBoundaryStencil(cell)
{
    boundaryPoints_.push_back(shape.nearestIntersect(cell_.centroid()));
    boundaryNormals_.push_back((boundaryPoints_.back() - cell_.centroid()).unitVec());

    auto hasNeighbourInSolid = [&shape](const Cell& cell) {
        for(const InteriorLink& nb: cell.neighbours())
            if(shape.isInside(nb.cell().centroid()))
                return true;

        return false;
    };

    for(const InteriorLink& nb: cell.neighbours())
        if(!shape.isInside(nb.cell().centroid()))
        {
            cellPoints_.push_back(nb.cell());

            if(hasNeighbourInSolid(nb.cell()))
            {
                boundaryPoints_.push_back(shape.nearestIntersect(nb.cell().centroid()));
                boundaryNormals_.push_back((boundaryPoints_.back() - cell_.centroid()).unitVec());
            }
        }

    for(const DiagonalCellLink& dg: cell.diagonals())
        if(!shape.isInside(dg.cell().centroid()))
        {
            cellPoints_.push_back(dg.cell());

            if(hasNeighbourInSolid(dg.cell()))
            {
                boundaryPoints_.push_back(shape.nearestIntersect(dg.cell().centroid()));
                boundaryNormals_.push_back((boundaryPoints_.back() - cell_.centroid()).unitVec());
            }
        }

    dMat_ = Matrix(6, cellPoints_.size() + boundaryPoints_.size());
    nMat_ = Matrix(6, cellPoints_.size() + boundaryPoints_.size());

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

    for(const Point2D& pt: boundaryPoints_)
    {
        dMat_(i, 0) = pt.x*pt.x;
        dMat_(i, 1) = pt.y*pt.y;
        dMat_(i, 2) = pt.x*pt.y;
        dMat_(i, 3) = pt.x;
        dMat_(i, 4) = pt.y;
        dMat_(i, 5) = 1.;

        Vector2D n = (boundaryPoints_.front() - cell_.centroid()).unitVec();

        nMat_(i, 0) = 2*n.x*pt.x;
        nMat_(i, 1) = 2*n.y*pt.y;
        nMat_(i, 2) = n.x*pt.y + n.y*pt.x;
        nMat_(i, 3) = n.x;
        nMat_(i, 4) = n.y;
        nMat_(i, 5) = 0.;

        i++;
    }

    dMat_ = dMat_*transpose(inverse(transpose(dMat_)*dMat_));
    nMat_ = nMat_*transpose(inverse(transpose(nMat_)*nMat_));
}
