#include "CellSearch.h"

void CellSearch::constructRTree()
{
    rTree_.clear();
    for(const Cell &cell: cellGroup_)
        rTree_.insert(Value(cell.centroid(), std::cref(cell)));
}

std::vector< Ref<const Cell> >  CellSearch::rangeSearch(const Circle &circle) const
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

std::vector< Ref<const Cell> > CellSearch::rangeSearch(const Polygon &pgn) const
{
    auto inPgn = [&pgn](const Value &val) { return pgn.isInside(val.first); };
    std::vector< Value > result;

    boost::geometry::model::box<Point2D> box;
    boost::geometry::envelope(pgn.boostPolygon(), box);

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(inPgn)
                 ,std::back_inserter(result));

    return getRefs(result);
}

std::vector< Ref<const Cell> > CellSearch::kNearestNeighbourSearch(const Point2D& point, size_t k) const
{
    std::vector< Value > result;

    rTree_.query(boost::geometry::index::nearest(point, k),
                 std::back_inserter(result));

    return getRefs(result);
}

//- Private methods

std::vector< Ref<const Cell> > CellSearch::getRefs(const std::vector<Value> &vals) const
{
    std::vector< Ref<const Cell> > refs;

    for(const Value &val: vals)
        refs.push_back(val.second);

    return refs;
}
