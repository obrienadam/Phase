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
    if(cellSet_.insert(std::make_pair(cell.id(), std::cref(cell))).second)
    {
        cells_.push_back(std::cref(cell));
        rTree_.insert(Value(cell.centroid(), cell.id()));
    }
}

void CellGroup::remove(const Cell &cell)
{
    auto entry = cellSet_.find(cell.id());

    if(entry != cellSet_.end())
    {
        cells_.erase(std::remove_if(cells_.begin(), cells_.end(), [&cell](const Cell &c){ return cell.id() == c.id(); }));
        cellSet_.erase(entry);
        rTree_.remove(Value(cell.centroid(), cell.id()));
    }
}

//- Searching
std::vector< Ref<const Cell> > CellGroup::rangeSearch(const Shape2D& shape) const
{
    auto insideShape = [&shape](const Value &val) { return shape.isInside(val.first); };
    std::vector<Value> result;

    boost::geometry::model::box<Point2D> box = shape.boundingBox();

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(insideShape),
                 std::back_inserter(result));

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
    auto it = cellSet_.find(cell.id());
    return it == cellSet_.end() ? false : true;
}

//- Protected

std::vector< Ref<const Cell> > CellGroup::getRefs(const std::vector< Value >& vals) const
{
    std::vector< Ref<const Cell> > refs;
    refs.reserve(vals.size());

    for(const Value &val: vals)
        refs.push_back(cellSet_.find(val.second)->second);

    return refs;
}
