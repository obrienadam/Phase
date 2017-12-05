#include "HighOrderStencil.h"

HighOrderStencil::HighOrderStencil(const Cell &cell, const ImmersedBoundary &ib)
{
    auto getIbObj = [&ib](const Cell& cell)->std::shared_ptr<const ImmersedBoundaryObject>
    {
        for(const CellLink& nb: cell.neighbours())
        {
            auto ibObj = ib.ibObj(nb.cell().centroid());
            if(ibObj)
                return ibObj;
        }

        return nullptr;
    };

    for(const CellLink& nb: cell.cellLinks())
    {
        if(ib.ibObj(nb.cell().centroid()))
            continue;

        auto ibObj = getIbObj(nb.cell());

        if(ibObj) //- is a compatibility/fluid point
        {
            compatCells_.push_back(nb.cell());
            compatPoints_.push_back(ibObj->nearestIntersect(nb.cell().centroid()));
        }
        else //- is just a fluid point
        {
            cells_.push_back(nb.cell());
        }

        //- This assumes that it is an IB point
        boundaryPoints_.push_back(ib.ibObj(cell.centroid())->nearestIntersect(cell.centroid()));
    }
}

HighOrderStencil::HighOrderStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
{
    auto isCompatCell = [&ibObj](const Cell& cell)
    {
        for(const CellLink& nb: cell.neighbours())
            if(ibObj.isInIb(nb.cell()))
                return true;
        return false;
    };

    for(const CellLink& nb: cell.cellLinks())
    {
        if(ibObj.isInIb(nb.cell()))
            continue;

        if(isCompatCell(nb.cell())) //- is a compatibility/fluid point
        {
            compatCells_.push_back(nb.cell());
            compatPoints_.push_back(ibObj.nearestIntersect(nb.cell().centroid()));
        }
        else //- is just a fluid point
        {
            cells_.push_back(nb.cell());
        }

        //- This assumes that it is an IB point
        boundaryPoints_.push_back(ibObj.nearestIntersect(cell.centroid()));
    }
}