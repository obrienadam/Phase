#ifndef INDEX_MAP_H
#define INDEX_MAP_H

#include "FiniteVolumeGrid2D.h"

class IndexMap
{
public:

    enum {INACTIVE = -1};

    IndexMap() {}

    IndexMap(const FiniteVolumeGrid2D& grid, Size nIndexSets = 1);

    void init(const FiniteVolumeGrid2D& grid, Size nIndexSets = 1);

    Index local(const Cell& cell, int set = 0) const;

    Index global(const Cell& cell, int set = 0) const;

    bool isLocallyActive(const Cell& cell) const;

    bool isGloballyActive(const Cell& cell) const;

private:

    Size nCells_ = 0, nIndexSets_ = 0;
    std::pair<Index, Index> range_ = std::make_pair(0, 0);
    std::vector<Index> localIndices_, globalIndices_;
};

#endif