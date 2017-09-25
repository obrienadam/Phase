#include <numeric>

#include "IndexMap.h"

IndexMap::IndexMap(const FiniteVolumeGrid2D &grid, Size nIndexSets)
{
    init(grid, nIndexSets);
}

void IndexMap::init(const FiniteVolumeGrid2D &grid, Size nIndexSets)
{
    nCells_ = grid.nCells();
    nIndexSets_ = nIndexSets;
    localIndices_.resize(nIndexSets*nCells_, INACTIVE);
    globalIndices_.resize(nIndexSets*nCells_, INACTIVE);

    auto nLocalUnknowns = grid.comm().allGather(nIndexSets_ * grid.nLocalActiveCells());
    range_.first = std::accumulate(nLocalUnknowns.begin(), nLocalUnknowns.begin() + grid.comm().rank(), 0);
    range_.second = range_.first + nLocalUnknowns[grid.comm().rank()];

    for (Size indexSet = 0; indexSet < nIndexSets_; ++indexSet)
    {
        Index localIndex = 0;
        for (const Cell &cell: grid.localActiveCells())
        {
            localIndices_[indexSet * nCells_ + cell.id()] = localIndex;
            globalIndices_[indexSet * nCells_ + cell.id()] = range_.first + localIndex++;
        }
    }

    //- Communicate global indices to other procs
    grid.sendMessages(globalIndices_, nIndexSets_);
}

Index IndexMap::local(const Cell &cell, int set) const
{
    return localIndices_[set * nCells_ + cell.id()];
}

Index IndexMap::global(const Cell &cell, int set) const
{
    return globalIndices_[set * nCells_ + cell.id()];
}

bool IndexMap::isLocallyActive(const Cell& cell) const
{
    return local(cell) != INACTIVE;
}

bool IndexMap::isGloballyActive(const Cell& cell) const
{
    return global(cell) != INACTIVE;
}
