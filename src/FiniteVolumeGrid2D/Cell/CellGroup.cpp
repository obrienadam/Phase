#include <algorithm>

#include "CellGroup.h"
#include "FiniteVolumeGrid2D.h"

void CellGroup::clear()
{
    cellSet_.clear();
    cells_.clear();
    rTree_.clear();
    geomRTree_.clear();
}

void CellGroup::push_back(Cell &cell)
{
    if (cellSet_.insert(std::make_pair(cell.id(), std::ref(cell))).second)
    {
        cells_.push_back(std::ref(cell));
        rTree_.insert(Value(cell.centroid(), cell.id()));
        geomRTree_.insert(std::make_pair(cell.shape().boundingBox(), cell.id()));
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
        geomRTree_.remove(std::make_pair(cell.shape().boundingBox(), cell.id()));
    }
}

void CellGroup::merge(CellGroup &other)
{
    for (Cell &cell: other.cells())
    {
        push_back(cell);
        other.remove(cell);
    }
}

//- Searching
std::vector<Ref<Cell> > CellGroup::cellCentersWithin(const Shape2D &shape)
{
    auto insideShape = [&shape](const Value &val) {
        return shape.isCovered(val.first);
    };

    std::vector<Value> result;

    boost::geometry::model::box<Point2D> box = shape.boundingBox();

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(insideShape),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<Cell> > CellGroup::cellCentersWithin(const Circle &circle)
{
    auto isInCircle = [&circle](const Value& val)
    {
        return circle.isInside(val.first);
    };

    std::vector<Value> result;
    rTree_.query(boost::geometry::index::covered_by(circle.boundingBox()) &&
                 boost::geometry::index::satisfies(isInCircle),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<Cell> > CellGroup::cellCentersWithin(const Box& box)
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::covered_by(box.boundingBox()),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<Cell> > CellGroup::cellNearestNeighbours(const Point2D &pt, size_t k)
{
    std::vector<Value> result;

    rTree_.query(boost::geometry::index::nearest(pt, k),
                 std::back_inserter(result));

    return getRefs(result);
}

std::vector<Ref<const Cell> > CellGroup::cellCentersWithin(const Shape2D &shape) const
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(shape.polygonize().boostPolygon()),
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

std::vector<Ref<Cell> > CellGroup::getRefs(const std::vector<Value> &vals)
{
    std::vector<Ref<Cell> > refs;
    refs.reserve(vals.size());

    for (const Value &val: vals)
        refs.push_back(cellSet_.find(val.second)->second);

    return refs;
}

std::vector<Ref<const Cell> > CellGroup::getRefs(const std::vector<Value> &vals) const
{
    std::vector<Ref<const Cell> > refs;
    refs.reserve(vals.size());

    for (const Value &val: vals)
        refs.push_back(std::cref(cellSet_.find(val.second)->second));

    return refs;
}
