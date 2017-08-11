#include "IndexMap.h"

IndexMap::IndexMap(const FiniteVolumeGrid2D &grid, Size nIndexSets)
    :
    nCells_(grid.nCells()),
    localIndices_(2 * nIndexSets, INACTIVE),
    globalIndices_(2 * nIndexSets, INACTIVE)
{
    const CellZone& localActiveCells = grid.localActiveCells();
    std::vector<Size> nIndicesProc = grid.comm().allGather(localActiveCells.size());

    for(Size &nIndices: nIndicesProc)
        nIndices *= nIndexSets;

    Size nIndicesThisProc = nIndicesProc[grid.comm().rank()];
    Size globalIndexStart = 0;

    for(int proc = 0; proc < grid.comm().rank(); ++proc)
        globalIndexStart += nIndicesProc[proc];

    for(Size indexSet = 0; indexSet < nIndexSets; ++indexSet)
    {
        Index idx = 0;
        for(const Cell& cell: localActiveCells)
            localIndices_[nCells_ * indexSet + cell.id()] = indexSet*nIndicesThisProc + idx++;

        idx = 0;

        for(const Cell& cell: localActiveCells)
            globalIndices_[nCells_ * indexSet + cell.id()] = globalIndexStart + indexSet*nIndicesThisProc + idx++;

        // grid.sendMessagesIt(globalIndices_.begin());
    }
}

Index IndexMap::local(const Cell& cell, int set) const
{
    return localIndices_[set*nCells_ + cell.id()];
}

Index IndexMap::global(const Cell& cell, int set) const
{
    return globalIndices_[set*nCells_ + cell.id()];
}
