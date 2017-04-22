#ifndef CELL_GROUP_H
#define CELL_GROUP_H

#include <vector>
#include <map>

#include <boost/geometry/index/rtree.hpp>

#include "Cell.h"
#include "Shape2D.h"
#include "Circle.h"
#include "Box.h"

typedef std::set<Ref<const Cell> > CellSet;

class CellGroup
{
public:

    CellGroup(const std::string &name, FiniteVolumeGrid2D &grid) : name_(name), grid_(grid)
    {}

    //- Cell group name
    const std::string &name() const
    { return name_; }

    const std::string &rename(const std::string &name)
    { return name_ = name; }

    //- Sizing
    virtual void clear();

    void reserve(size_t size)
    { cells_.reserve(size); }

    size_t size() const
    { return cells_.size(); }

    bool empty() const
    { return cells_.empty(); }

    //- Data access
    std::vector<Ref<Cell> > cells() const
    { return cells_; }

    //- Adding/removing cells
    virtual void push_back(Cell &cell);

    virtual void remove(const Cell &cell);

    virtual void merge(CellGroup &other);

    //- Searching
    std::vector<Ref<Cell> > cellCentersWithin(const Shape2D &shape);

    std::vector<Ref<Cell> > cellCentersWithin(const Circle &circle);

    std::vector<Ref<Cell> > cellCentersWithin(const Box &box);

    std::vector<Ref<Cell> > cellsOverlapping(const Shape2D &shape);

    std::vector<Ref<Cell> > cellNearestNeighbours(const Point2D &pt, size_t k);

    std::vector<Ref<const Cell> > cellCentersWithin(const Shape2D &shape) const;

    std::vector<Ref<const Cell> > cellCentersWithin(const Circle &circle) const;

    std::vector<Ref<const Cell> > cellCentersWithin(const Box &box) const;

    std::vector<Ref<const Cell> > cellsOverlapping(const Shape2D &shape) const;

    std::vector<Ref<const Cell> > cellNearestNeighbours(const Point2D &pt, size_t k) const;

    //- Iterators
    std::vector<Ref<Cell> >::iterator begin()
    { return cells_.begin(); }

    std::vector<Ref<Cell> >::iterator end()
    { return cells_.end(); }

    std::vector<Ref<Cell> >::const_iterator begin() const
    { return cells_.begin(); }

    std::vector<Ref<Cell> >::const_iterator end() const
    { return cells_.end(); }

    bool isInGroup(const Cell &cell) const;

protected:

    typedef std::pair<Point2D, Label> Value;
    typedef boost::geometry::index::rtree<Value, boost::geometry::index::quadratic<16> > Rtree;

    std::vector<Ref<Cell> > getRefs(const std::vector<Value> &vals);

    std::vector<Ref<const Cell> > getRefs(const std::vector<Value> &vals) const;

    std::string name_;

    FiniteVolumeGrid2D &grid_;

    std::map<Label, Ref<Cell> > cellSet_; // Allows cell lookup via an id
    std::vector<Ref<Cell> > cells_; // Used for faster iteration over all cells

    std::map<Label, Ref<Node> > nodeSet_;
    std::vector<Ref<Node> > nodes_;

    Rtree rTree_;
    boost::geometry::index::rtree<std::pair<boost::geometry::model::box<Point2D>, Label>, boost::geometry::index::quadratic<16>> geomRTree_;
};

#endif
