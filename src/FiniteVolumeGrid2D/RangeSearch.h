#ifndef RANGE_SEARCH_H
#define RANGE_SEARCH_H

#include "Geometry.h"
#include "FiniteVolumeGrid2D.h"

class RangeSearch
{
public:

    RangeSearch() {}
    RangeSearch(const FiniteVolumeGrid2D& grid, Scalar radius);

    void search(const FiniteVolumeGrid2D& grid, Scalar radius);

    const std::vector< Ref<const Cell> >& getResult(size_t id) const { return result_[id]; }

private:

    std::vector< std::vector< Ref<const Cell> > > result_;

};

#endif
