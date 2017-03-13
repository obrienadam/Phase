#ifndef SEARCH_H
#define SEARCH_H

#include <boost/geometry/index/rtree.hpp>

#include "Circle.h"

class Search
{
public:

    Search()
    {}

    Search(const std::vector<Point2D> &points)
    { constructRTree(points); }

    void constructRTree(const std::vector<Point2D> &points);

    void add(const Point2D &point, Label label);

    void clear()
    { rTree_.clear(); }

    std::vector<Label> rangeSearch(const Circle &circle) const;

    std::vector<Label> kNearestNeighbourSearch(const Point2D &point, size_t k) const;

    std::vector<Label> rangeSearch(const Shape2D &shape) const;

private:

    boost::geometry::index::rtree<std::pair<Point2D, Label>, boost::geometry::index::quadratic<32> > rTree_;
};

#endif
