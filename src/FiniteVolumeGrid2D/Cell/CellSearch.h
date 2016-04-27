#ifndef CELL_SEARCH_H
#define CELL_SEARCH_H

#include "Geometry.h"
#include "FiniteVolumeGrid2D.h"
#include "Circle.h"

std::vector< Ref<const Cell> > rangeSearch(const FiniteVolumeGrid2D &grid, const Circle& circle);
std::vector< std::vector< Ref<const Cell> > > rangeSearch(const FiniteVolumeGrid2D& grid, Scalar radius);

std::vector< Ref<const Cell> > kNearestNeighbourSearch(const FiniteVolumeGrid2D &grid, const Cell& cell, size_t k);
std::vector< std::vector< Ref<const Cell> > > kNearestNeighbourSearch(const FiniteVolumeGrid2D& grid, size_t k);

#endif
