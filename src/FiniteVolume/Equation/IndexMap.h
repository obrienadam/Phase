#ifndef INDEX_MAP_H
#define INDEX_MAP_H

#include "FiniteVolumeGrid2D.h"

class IndexMap
{
public:

    enum
    {
        INACTIVE = -1
    };

    IndexMap()
    {}

    IndexMap(const FiniteVolumeGrid2D &grid, Size nIndices = 1)
    { init(grid, nIndices); }

    void init(const FiniteVolumeGrid2D &grid, Size nIndices);

    void update(const FiniteVolumeGrid2D &grid);

    Index local(const Cell &cell, Label indexNo = 0) const
    { return localIndices_[indexNo * nCells_ + cell.id()]; }

    Index global(const Cell &cell, Label indexNo = 0) const
    { return globalIndices_[indexNo * nCells_ + cell.id()]; }

    bool isActive(const Cell &cell) const
    { return globalIndices_[cell.id()] != -1; }

    const std::pair<Index, Index> &ownershipRange() const
    { return ownershipRange_; }

    Index minGlobalIndex() const
    { return ownershipRange_.first; }

    Index maxGlobalIndex() const
    { return ownershipRange_.second - 1; }

private:

    Size nCells_, nIndices_;
    std::pair<Index, Index> ownershipRange_;
    std::vector<Index> localIndices_, globalIndices_;
};

#endif