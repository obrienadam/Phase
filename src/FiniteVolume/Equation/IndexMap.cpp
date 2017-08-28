#include "IndexMap.h"

IndexMap::IndexMap(const FiniteVolumeGrid2D &grid, Size nIndexSets)
    :
    nCells_(grid.nCells()),
    nIndexSets_(nIndexSets),
    localIndices_(nIndexSets_ * grid.nCells(), INACTIVE),
    globalIndices_(nIndexSets_ * grid.nCells(), INACTIVE)
{
    auto nUnknowns = grid.comm().allGather(nIndexSets_ * grid.nLocalActiveCells());

    for(int proc = 0; proc < grid.comm().rank(); ++proc)
        range_.first += nUnknowns[proc];
    range_.second = range_.first + nUnknowns[grid.comm().rank()];

    for(Size indexSet = 0; indexSet < nIndexSets_; ++indexSet)
    {
        Index idx = 0;
        for(const Cell& cell: grid.localActiveCells())
        {
            localIndices_[indexSet * nCells_ + cell.id()] = idx;
            globalIndices_[indexSet * nCells_ + cell.id()] = range_.first + idx++;
        }
    }

    //- Communicate global indices to other procs
    grid.sendMessages(globalIndices_, nIndexSets_);
}

Index IndexMap::local(const Cell& cell, int set) const
{
    return localIndices_[set*nCells_ + cell.id()];
}

Index IndexMap::global(const Cell& cell, int set) const
{
    return globalIndices_[set*nCells_ + cell.id()];
}
