#include "DirectForcingImmersedBoundaryLeastSquaresQuadraticStencil.h"

DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::LeastSquaresQuadraticStencil(const Cell &cell,
                                                                                          const DirectForcingImmersedBoundary &ib)
{
    _cells.reserve(8);
    _compatPts.reserve(8);

    for(const CellLink &nb: cell.cellLinks())
    {
        if(!ib.globalSolidCells().isInSet(nb.cell()))
        {
            _cells.push_back(&nb.cell());

            if(ib.globalIbCells().isInSet(nb.cell()))
                _compatPts.push_back(CompatPoint(nb.cell(), *ib.nearestIbObj(nb.cell().centroid())));
        }
    }

    _compatPts.push_back(CompatPoint(cell, *ib.nearestIbObj(cell.centroid())));
}
