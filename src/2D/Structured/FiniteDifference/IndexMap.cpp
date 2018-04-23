#include "IndexMap.h"

IndexMap::IndexMap(const StructuredGrid2D &grid, int nIndexSets)
{
    _nNodes = grid.nNodes();
    _nIndexSets = nIndexSets;

    _localIndices.resize(_nNodes * _nIndexSets);
    std::iota(_localIndices.begin(), _localIndices.end(), 0);

    _globalIndices.resize(_nNodes * _nIndexSets);
    std::iota(_globalIndices.begin(), _globalIndices.end(), 0);
}