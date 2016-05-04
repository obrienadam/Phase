#ifndef CELL_SEARCH_H
#define CELL_SEARCH_H

#include <boost/geometry/index/rtree.hpp>

#include "Circle.h"
#include "CellGroup.h"

class CellSearch
{
public:

    CellSearch(const CellGroup &cellGroup) : cellGroup_(cellGroup) { constructRTree(); }

    void constructRTree();

    std::vector< Ref<const Cell> > rangeSearch(const Circle &circle) const;
    std::vector< Ref<const Cell> > rangeSearch(const Polygon &pgn) const;
    std::vector< Ref<const Cell> > kNearestNeighbourSearch(const Point2D& point, size_t k) const;

    const CellGroup& cellGroup() const { return cellGroup_; }

private:

    typedef std::pair< Point2D, Ref< const Cell > > Value;

    std::vector< Ref<const Cell> > getRefs(const std::vector< Value >& vals) const;

    boost::geometry::index::rtree< Value, boost::geometry::index::quadratic<32> > rTree_;

    const CellGroup &cellGroup_;

};

#endif
