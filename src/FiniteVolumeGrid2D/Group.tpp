#include "Group.h"

template<class T>
void Group<T>::clear()
{
    items_.clear();
    itemSet_.clear();
    rTree_.clear();
}

template<class T>
void Group<T>::add(const T &item)
{
    if (itemSet_.insert(std::make_pair(item.id(), std::cref(item))).second)
    {
        items_.push_back(std::cref(item));
        rTree_.insert(Value(item.centroid(), item.id()));
    }
}

template<class T>
void Group<T>::add(const Group<T> &other)
{
    for (const T &item: other)
        add(item);
}

template<class T>
void Group<T>::remove(const T &item)
{
    auto entry = itemSet_.find(item.id());

    if (entry != itemSet_.end())
    {
        itemSet_.erase(entry);
        items_.erase(std::remove_if(items_.begin(), items_.end(), [&item](const T &i) { return item.id() == i.id(); }));
        rTree_.remove(Value(item.centroid(), item.id()));
    }
}

template<class T>
void Group<T>::remove(const Group<T> &other)
{
    for (const T &item: other)
        remove(item);
}

template<class T>
Group<T> &Group<T>::operator+=(const Group<T> &rhs)
{
    if (this != &rhs)
        add(rhs.begin(), rhs.end());
    return *this;
}

template<class T>
Group<T> &Group<T>::operator-=(const Group<T> &rhs)
{
    remove(rhs);
    return *this;
}

template<class T>
Group<T> Group<T>::intersection(const Group<T> &other) const
{
    Group<T> result;

    for (const T &item: items_)
        if (other.isInGroup(item))
            result.add(item);

    return result;
}

template<class T>
Group<T> Group<T>::difference(const Group<T> &other) const
{
    Group<T> result;

    for (const T &item: items_)
        if (!other.isInGroup(item))
            result.add(item);

    return result;
}

template<class T>
Group<T> Group<T>::symmetricDifference(const Group<T> &other) const
{
    Group<T> result;

    for (const T &item: items_)
        if (!other.isInGroup(item))
            result.add(item);

    for (const T &item: other.items_)
        if (!isInGroup(item))
            result.add(item);

    return result;
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsWithin(const Shape2D &shape) const
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(shape.polygonize().boostRing()),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsWithin(const Circle &circle) const
{
    auto isInCircle = [&circle](const Value &val) {
        return circle.isInside(val.first);
    };

    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(circle.boundingBox()) &&
                 boost::geometry::index::satisfies(isInCircle),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsWithin(const Box &box) const
{
    std::vector<Value> result;
    rTree_.query(boost::geometry::index::within(box.boundingBox()),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsCoveredBy(const Circle &circle) const
{
    std::vector<Value> result;

    auto isCoveredByCircle = [&circle](const Value &val) {
        return circle.isCovered(val.first);
    };

    rTree_.query(boost::geometry::index::covered_by(circle.boundingBox()) &&
                 boost::geometry::index::satisfies(isCoveredByCircle),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsCoveredBy(const Box &box) const
{
    std::vector<Value> result;

    rTree_.query(boost::geometry::index::covered_by(box.boundingBox()),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::nearestItems(const Point2D &pt, size_t k) const
{
    std::vector<Value> result;

    rTree_.query(boost::geometry::index::nearest(pt, k),
                 std::back_inserter(result));

    return getRefs(result);
}

template<class T>
const T &Group<T>::nearestItem(const Point2D &pt) const
{
    return nearestItems(pt, 1)[0];
}

template<class T>
bool Group<T>::isInGroup(const T &item) const
{
    return itemSet_.find(item.id()) != itemSet_.end();
}

//- Private helper methods
template<class T>
std::vector<Ref<const T> > Group<T>::getRefs(const std::vector<Group::Value> &vals) const
{
    std::vector<Ref<const T> > refs;
    refs.reserve(vals.size());

    std::transform(vals.begin(), vals.end(), std::back_inserter(refs),
                   [this](const Value &v) { return itemSet_.find(v.second)->second; });

    return refs;
}

//- Some useful operators
template<class T>
Group<T> operator+(Group<T> lhs, const Group<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Group<T> operator-(Group<T> lhs, const Group<T> &rhs)
{
    lhs -= rhs;
    return lhs;
}