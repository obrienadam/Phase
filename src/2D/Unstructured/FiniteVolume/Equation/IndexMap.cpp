#include <numeric>

#include "IndexMap.h"

void IndexMap::init(const FiniteVolumeGrid2D &grid, Size nIndices)
{
    nCells_ = grid.nCells();
    nIndices_ = nIndices;

    localIndices_.resize(nIndices_ * nCells_);
    globalIndices_.resize(nIndices_ * nCells_);

    update(grid);
}

void IndexMap::update(const FiniteVolumeGrid2D &grid)
{
    std::fill(localIndices_.begin(), localIndices_.end(), INACTIVE);
    std::fill(globalIndices_.begin(), globalIndices_.end(), INACTIVE);

    std::vector<Size> nLocalActiveCells = grid.comm().allGather(grid.localCells().size());

    ownershipRange_.first =
            nIndices_ * std::accumulate(nLocalActiveCells.begin(), nLocalActiveCells.begin() + grid.comm().rank(), 0);
    ownershipRange_.second = ownershipRange_.first + nIndices_ * nLocalActiveCells[grid.comm().rank()];

    Index localIndex = 0;

    for (Size indexNo = 0; indexNo < nIndices_; ++indexNo)
        for (const Cell &cell: grid.localCells())
        {
            localIndices_[indexNo * nCells_ + cell.id()] = localIndex;
            globalIndices_[indexNo * nCells_ + cell.id()] = ownershipRange_.first + localIndex++;
        }

    //- Communicate global indices to other procs
    grid.sendMessages(globalIndices_, nIndices_);
}
