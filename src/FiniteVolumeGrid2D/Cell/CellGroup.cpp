#include <algorithm>

#include "CellGroup.h"
#include "FiniteVolumeGrid2D.h"

void CellGroup::clear()
{
    cellSet_.clear();
    cells_.clear();
    rTree_.clear();
}

void CellGroup::push_back(const Cell &cell)
{
    if (cellSet_.insert(std::make_pair(cell.id(), std::ref(cell))).second)
    {
        cells_.push_back(std::cref(cell));
        rTree_.insert(Value(cell.centroid(), cell.id()));
    }
}

void CellGroup::remove(const Cell &cell)
{
    auto entry = cellSet_.find(cell.id());

    if (entry != cellSet_.end())
    {
        cells_.erase(
                std::remove_if(cells_.begin(), cells_.end(), [&cell](const Cell &c) { return cell.id() == c.id(); }));
        cellSet_.erase(entry);
        rTree_.remove(Value(cell.centroid(), cell.id()));
    }
}

void CellGroup::merge(CellGroup &other)
{
    for (const Cell &cell: other.cells())
    {
        push_back(cell);
        other.remove(cell);
    }
}

//- Searching
std::vector<Ref<const Cell> > CellGroup::cellCentersWithin(const Shape2D &shape) const
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(shape.polygonize().boostRing()),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<const Cell> > CellGroup::cellCentersWithin(const Circle &circle) const
{
    auto isInCircle = [&circle](const Value& val)
    {
        return circle.isInside(val.first);
    };

    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(circle.boundingBox()) &&
                 boost::geometry::index::satisfies(isInCircle),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<const Cell> > CellGroup::cellCentersWithin(const Box &box) const
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(box.boundingBox()),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<const Cell> > CellGroup::cellNearestNeighbours(const Point2D &pt, size_t k) const
{
    std::vector<Value> result;

    rTree_.query(boost::geometry::index::nearest(pt, k),
                 std::back_inserter(result));

    return getRefs(result);
}

bool CellGroup::isInGroup(const Cell &cell) const
{
    return cellSet_.find(cell.id()) != cellSet_.end();
}

//- Protected

std::vector<Ref<const Cell> > CellGroup::getRefs(const std::vector<Value> &vals) const
{
    std::vector<Ref<const Cell> > refs;
    refs.reserve(vals.size());

    for (const Value &val: vals)
        refs.push_back(std::cref(cellSet_.find(val.second)->second));

    return refs;
}
