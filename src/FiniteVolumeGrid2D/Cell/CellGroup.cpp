#include <algorithm>

#include "CellGroup.h"

void CellGroup::clear()
{
    cellSet_.clear();
    cells_.clear();
    rTree_.clear();
}

void CellGroup::push_back(const Cell &cell)
{
    if(cellSet_.insert(std::make_pair(&cell, cells_.size())).second)
    {
        cells_.push_back(std::cref(cell));
        rTree_.insert(Value(cell.centroid(), std::cref(cell)));
    }
}

void CellGroup::remove(const Cell &cell)
{
    auto entry = cellSet_.find(&cell);

    if(entry != cellSet_.end())
    {
        cells_.erase(cells_.begin() + entry->second);
        cellSet_.erase(entry);

        size_t idx = 0;
        for(const Cell &cell: cells_) // Must reset remaining cell indices, very expensive!
            cellSet_[&cell] = idx++;

        rTree_.clear();
        for(const Cell &cell: cells())
            rTree_.insert(Value(cell.centroid(), std::cref(cell)));
    }
}

//- Searching
std::vector< Ref<const Cell> > CellGroup::rangeSearch(const Circle& circle) const
{
    auto inCircle = [&circle](const Value &val) { return circle.isInside(val.first); };
    std::vector< Value > result;

    boost::geometry::model::box<Point2D> box(Point2D(circle.centroid().x - circle.radius(),
                                                     circle.centroid().y - circle.radius()),
                                             Point2D(circle.centroid().x + circle.radius(),
                                                     circle.centroid().y + circle.radius()));

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(inCircle)
                 ,std::back_inserter(result));

    return getRefs(result);
}

std::vector< Ref<const Cell> > CellGroup::kNearestNeighbourSearch(const Point2D& pt, size_t k) const
{
    std::vector< Value > result;

    rTree_.query(boost::geometry::index::nearest(pt, k),
                 std::back_inserter(result));

    return getRefs(result);
}

bool CellGroup::isInGroup(const Cell &cell) const
{
    auto it = cellSet_.find(&cell);

    return it == cellSet_.end() ? false : true;
}

//- Protected

std::vector< Ref<const Cell> > CellGroup::getRefs(const std::vector< Value >& vals) const
{
    std::vector< Ref<const Cell> > refs;
    refs.reserve(vals.size());

    for(const Value &val: vals)
        refs.push_back(val.second);

    return refs;
}
