#ifndef PHASE_INDEX_MAP_H
#define PHASE_INDEX_MAP_H

#include "StructuredGrid2D/StructuredGrid2D.h"

class IndexMap
{
public:

    enum
    {
        INACTIVE = -1
    };

    IndexMap(const StructuredGrid2D &grid, int nIndexSets);

    Size size() const
    { return _nNodes * _nIndexSets; }

    Index local(const StructuredGrid2D::Node &node, int setNo) const
    { return _localIndices[setNo * _nNodes + node.id()]; }

    Index global(const StructuredGrid2D::Node &node, int setNo) const
    { return _globalIndices[setNo * _nNodes + node.id()]; }

protected:

    Size _nNodes, _nIndexSets;

    std::vector<Index> _localIndices, _globalIndices;

};


#endif //PHASE_INDEXMAP_H
