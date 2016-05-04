#ifndef SEARCH_H
#define SEARCH_H

#include <boost/geometry/index/rtree.hpp>

#include "Circle.h"

template <class T>
class Search
{
public:

    Search(std::vector<T>& items) : items_(items) {}

    void constructRTree();

    std::vector< Ref<const T> > rangeSearch(const Circle &circle) const;
    std::vector< Ref<const T> > kNearestNeighbourSearch(const Point2D& point, size_t k) const;

private:

    typedef std::pair< Point2D, Ref< const T > > Value;

    std::vector< Ref<const T> > getRefs(const std::vector< Value >& vals) const;

    boost::geometry::index::rtree< Value, boost::geometry::index::quadratic<32> > rTree_;
    const std::vector<T>& items_;

};

#include "Search.tpp"

#endif
