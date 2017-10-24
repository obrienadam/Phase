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

    namespace bgi = boost::geometry::index;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
        {
            const Circle &c = static_cast<const Circle &>(shape);
            auto covered = [&c](const Value &v) { return c.isInside(v.first); };
            rTree_.query(bgi::within(c.boundingBox()) && bgi::satisfies(covered), std::back_inserter(result));
        }
            break;

        case Shape2D::BOX:
            rTree_.query(bgi::within(shape.boundingBox()), std::back_inserter(result));
            break;
        case Shape2D::POLYGON:
            rTree_.query(bgi::within(static_cast<const Polygon &>(shape).boostRing()), std::back_inserter(result));
            break;
    }


    return getRefs(result);
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsCoveredBy(const Shape2D &shape) const
{
    std::vector<Value> result;

    namespace bgi = boost::geometry::index;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
        {
            const Circle &c = static_cast<const Circle &>(shape);
            auto covered = [&c](const Value &v) { return c.isCovered(v.first); };
            rTree_.query(bgi::covered_by(c.boundingBox()) && bgi::satisfies(covered), std::back_inserter(result));
        }
            break;

        case Shape2D::BOX:
            rTree_.query(bgi::covered_by(shape.boundingBox()), std::back_inserter(result));
            break;
        case Shape2D::POLYGON:
            rTree_.query(bgi::covered_by(static_cast<const Polygon &>(shape).boostRing()), std::back_inserter(result));
            break;
    }

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
std::vector<Ref<const T>> Group<T>::nearestItems(const Shape2D &shape, size_t k) const
{
    std::vector<Value> result;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:

            rTree_.query(boost::geometry::index::nearest(shape.centroid(), k),
                         std::back_inserter(result));

            break;

        default:
            throw Exception("Group<T>", "nearestItems", "shape type not supported.");
    }

    return getRefs(result);
}

template<class T>
const T &Group<T>::nearestItem(const Point2D &pt) const
{
    return nearestItems(pt, 1)[0];
}

template<class T>
const T &Group<T>::nearestItem(const Shape2D &shape) const
{
    return nearestItems(shape, 1)[0];
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