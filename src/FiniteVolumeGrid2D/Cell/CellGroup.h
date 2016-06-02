#ifndef CELL_GROUP_H
#define CELL_GROUP_H

#include <vector>
#include <map>

#include <boost/geometry/index/rtree.hpp>

#include "Cell.h"
#include "Circle.h"

typedef std::set< Ref<const Cell> > CellSet;

class CellGroup
{
public:

    CellGroup(const std::string& name = "N/A") : name_(name) {}

    //- Cell group name
    const std::string& name() const { return name_; }
    const std::string& rename(const std::string& name) { return name_ = name; }

    //- Sizing
    virtual void clear();
    void reserve(size_t size) { cells_.reserve(size); }
    size_t size() const { return cells_.size(); }

    //- Data access
    const std::vector< Ref<const Cell> >& cells() const { return cells_; }

    //- Adding/removing cells
    virtual void push_back(const Cell &cell);
    virtual void remove(const Cell &cell);

    //- Searching
    std::vector< Ref<const Cell> > rangeSearch(const Circle& circle) const;
    std::vector< Ref<const Cell> > kNearestNeighbourSearch(const Point2D& pt, size_t k) const;

    //- Iterators
    std::vector< Ref<const Cell> >::iterator begin() { return cells_.begin(); }
    std::vector< Ref<const Cell> >::iterator end() { return cells_.end(); }
    std::vector< Ref<const Cell> >::const_iterator begin() const { return cells_.begin(); }
    std::vector< Ref<const Cell> >::const_iterator end() const { return cells_.end(); }

    bool isInGroup(const Cell &cell) const;

protected:

    typedef std::pair< Point2D, Ref< const Cell > > Value;

    std::vector< Ref<const Cell> > getRefs(const std::vector< Value >& vals) const;

    std::string name_;

    std::map< const Cell*, size_t> cellSet_;
    std::vector< Ref<const Cell> > cells_;

    boost::geometry::index::rtree< Value, boost::geometry::index::quadratic<32> > rTree_;
};

#endif
