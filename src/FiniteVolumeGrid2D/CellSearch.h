#ifndef CELL_SEARCH_H
#define CELL_SEARCH_H

#include "Geometry.h"
#include "FiniteVolumeGrid2D.h"

std::vector< std::vector< Ref<const Cell> > > rangeSearch(const FiniteVolumeGrid2D& grid, Scalar radius);
std::vector< std::vector< Ref<const Cell> > > kNearestNeighbourSearch(const FiniteVolumeGrid2D& grid, size_t k);

#endif
