#include "ForcingCellStencil.h"
#include "Exception.h"

ForcingCellStencil::ForcingCellStencil(const Cell &cell, const Shape2D &shape, const CellGroup &cellGroup)
    :
      ImmersedBoundaryStencil(cell)
{
    xc_ = shape.nearestIntersect(cell_.centroid());
    auto cells = cellGroup.cellNearestNeighbours(cell_.centroid(), 4);

    //- Remove the stencil cell from the nearest neighbour result
    cells.erase(std::remove_if(cells.begin(), cells.end(), [&cell](const Cell& nb){
        return nb.id() == cell.id();
    }), cells.end());

    std::vector<std::vector<Ref<const Cell>>> candidates;

    candidates.push_back({cells[0], cells[1]});
    candidates.push_back({cells[1], cells[2]});
    candidates.push_back({cells[0], cells[2]});

    iCells_ = candidates.front();
    Polygon bestStencil = {xc_, iCells_[0].get().centroid(), iCells_[1].get().centroid()};

    for(const auto& candidate: candidates)
    {
        Polygon stencil = {xc_, candidate[0].get().centroid(), candidate[1].get().centroid()};

        if(stencil.isInside(cell.centroid()) && stencil.perimeter() < bestStencil.perimeter())
        {
            iCells_ = candidate;
            bestStencil = stencil;
        }
    }

    interpolator_ = LinearInterpolation({xc_, iCells_[0].get().centroid(), iCells_[1].get().centroid()});
    std::vector<Scalar> coeffs = interpolator_(cell.centroid());

    bCoeff_ = coeffs[0];
    iCoeffs_ = {coeffs[1], coeffs[2]};

    coeffs = interpolator_.derivative(cell.centroid(), (xc_ - cell.centroid()).unitVec(), {iCells_[0].get().centroid(), iCells_[1].get().centroid()});

    bnCoeff_ = coeffs[0];
    nCoeffs_ = {coeffs[1], coeffs[2]};
}
