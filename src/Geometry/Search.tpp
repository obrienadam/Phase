#include "Search.h"

template<class T>
void Search<T>::constructRTree()
{
    rTree_.clear();
    for(const T& item: items_)
        rTree_.insert(Value(Point2D(item), std::cref(item)));
}

template<class T>
std::vector< Ref<const T> > Search<T>::rangeSearch(const Circle &circle) const
{
    std::vector< Value > result;

    auto inCircle = [&circle](const Value &val) { return circle.isInside(val.first); };

    boost::geometry::model::box<Point2D> box(Point2D(circle.centroid().x - circle.radius(),
                                                     circle.centroid().y - circle.radius()),
                                             Point2D(circle.centroid().x + circle.radius(),
                                                     circle.centroid().y + circle.radius()));

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(inCircle),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector< Ref<const T> > Search<T>::kNearestNeighbourSearch(const Point2D& point, size_t k) const
{
    std::vector< Value > result;

    rTree_.query(boost::geometry::index::nearest(point, k),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
void Search<T>::clear()
{
    rTree_.clear();
}

//- Private methods

template<class T>
std::vector< Ref<const T> > Search<T>::getRefs(const std::vector< Value >& vals) const
{
    std::vector< Ref<const T> > refs;

    for(const Value &val: vals)
        refs.push_back(val.second);

    return refs;
}
