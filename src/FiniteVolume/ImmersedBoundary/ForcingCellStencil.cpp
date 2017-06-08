#include "ForcingCellStencil.h"
#include "Exception.h"

ForcingCellStencil::ForcingCellStencil(const Cell &cell, const Shape2D &shape, const CellGroup &cellGroup)
    :
      ImmersedBoundaryStencil(cell)
{
    xc_ = shape.nearestIntersect(cell_.centroid());
    auto cells = cellGroup.nearestItems(cell_.centroid(), 3);

    Polygon stencils[3] = {
        Polygon({xc_, cells[0].get().centroid(), cells[1].get().centroid()}),
        Polygon({xc_, cells[0].get().centroid(), cells[2].get().centroid()}),
        Polygon({xc_, cells[1].get().centroid(), cells[2].get().centroid()})
    };

    Scalar minPerimeter = std::numeric_limits<Scalar>::infinity();
    int best = -1;

    for(int i = 0; i < 3; ++i)
    {
        Scalar perimeter = stencils[i].perimeter();

        if(stencils[i].isInside(cell.centroid()) && perimeter < minPerimeter)
        {
            best = i;
            minPerimeter = stencils[i].perimeter();
        }
    }

    switch(best)
    {
    case 0:
        nbCells_ = {cells[0], cells[1]};
        break;
    case 1:
        nbCells_ = {cells[0], cells[2]};
        break;
    case 2:
        nbCells_ = {cells[1], cells[2]};
        break;
    default:
        nbCells_ = {cells[0], cells[1]};;
    }

    Matrix dMat(3, 3);
    Matrix nMat(3, 3);

    Vector2D n = (cell.centroid() - xc_).unitVec();

    dMat(0, 0) = 0.;
    dMat(0, 1) = 0.;
    dMat(0, 2) = 1.;

    nMat(0, 0) = n.x;
    nMat(0, 1) = n.y;
    nMat(0, 2) = 0.;

    for(int i = 1; i < 3; ++i)
    {
        Point2D pt = nbCells_[i - 1].get().centroid() - xc_;

        dMat(i, 0) = pt.x;
        dMat(i, 1) = pt.y;
        dMat(i, 2) = 1.;

        nMat(i, 0) = pt.x;
        nMat(i, 1) = pt.y;
        nMat(i, 2) = 1.;
    }

    dMat.invert().transpose();
    nMat.invert().transpose();

    Matrix x(3, 1);

    Point2D pt = cell.centroid() - xc_;

    x = {
        pt.x,
        pt.y,
        1.
    };

    auto coeffs = dMat*x;

    dirichletBoundaryCoeff_ = coeffs(0, 0);
    dirichletCellCoeffs_ = {coeffs(1, 0), coeffs(2, 0)};

    coeffs = nMat*x;

    neumannBoundaryCoeff_ = coeffs(0, 0);
    neumannCellCoeffs_ = {coeffs(1, 0), coeffs(2, 0)};
}
