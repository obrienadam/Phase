#ifndef INDEX_MAP_H
#define INDEX_MAP_H

#include "FiniteVolumeGrid2D.h"

class IndexMap
{
public:

    enum {INACTIVE = -1};

    IndexMap(const FiniteVolumeGrid2D& grid, Size nIndexSets);

    Index local(const Cell& cell, int set = 0) const;

    Index global(const Cell& cell, int set = 0) const;

private:

    Size nCells_ = 0;
    std::vector<Index> localIndices_, globalIndices_;
};

#endif