#include "Search.h"

void Search::constructRTree(const std::vector<Point2D> &points)
{
    std::vector<std::pair<Point2D, Label>> values(points.size());
    Label id = 0;
    std::transform(points.begin(), points.end(),
                   values.begin(),
                   [&id](const Point2D &point) { return std::make_pair(point, id++); });

    rTree_.clear();
    rTree_.insert(values.begin(), values.end());
}

void Search::add(const Point2D &point, Label label)
{
    rTree_.insert(std::make_pair(point, label));
}

std::vector<Label> Search::rangeSearch(const Circle &circle) const
{
    std::vector<std::pair<Point2D, Label>> result;

    auto inCircle = [&circle](const std::pair<Point2D, Label> &val) { return circle.isInside(val.first); };

    boost::geometry::model::box<Point2D> box(Point2D(circle.centroid().x - circle.radius(),
                                                     circle.centroid().y - circle.radius()),
                                             Point2D(circle.centroid().x + circle.radius(),
                                                     circle.centroid().y + circle.radius()));

    rTree_.query(boost::geometry::index::within(box)
                 && boost::geometry::index::satisfies(inCircle),
                 std::back_inserter(result));

    std::vector<Label> ids(result.size());
    std::transform(result.begin(), result.end(), ids.begin(),
                   [](const std::pair<Point2D, Label> &value) { return value.second; });

    return ids;
}

std::vector<Label> Search::kNearestNeighbourSearch(const Point2D &point, size_t k) const
{
    std::vector<std::pair<Point2D, Label>> result;

    rTree_.query(boost::geometry::index::nearest(point, k),
                 std::back_inserter(result));

    std::vector<Label> ids(result.size());
    std::transform(result.begin(), result.end(), ids.begin(),
                   [](const std::pair<Point2D, Label> &value) { return value.second; });

    return ids;
}

std::vector<Label> Search::rangeSearch(const Shape2D &shape) const
{
    auto insideShape = [&shape](const std::pair<Point2D, Label> &val) {
        return shape.isCovered(val.first);
    };


    std::vector<std::pair<Point2D, Label>> result;
    rTree_.query(boost::geometry::index::within(shape.boundingBox())
                 && boost::geometry::index::satisfies(insideShape),
                 std::back_inserter(result));

    std::vector<Label> ids(result.size());
    std::transform(result.begin(), result.end(), ids.begin(),
                   [](const std::pair<Point2D, Label> &value) { return value.second; });

    return ids;
}
