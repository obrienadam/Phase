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
    if (itemSet_.insert(std::cref(item)).second)
    {
        items_.push_back(std::cref(item));
        rTree_.insert(std::cref(item));
    }
}

template<class T>
void Group<T>::add(const Group<T> &items)
{
    items_.reserve(items_.size() + items.size());
    itemSet_.reserve(itemSet_.size() + items.size());

    for (const T &item: items)
        add(item);
}

template<class T>
void Group<T>::remove(const T &item)
{
    if (itemSet_.erase(std::cref(item)))
    {
        items_.erase(
                std::find_if(items_.begin(), items_.end(), [&item](const T &arg) { return item.id() == arg.id(); })
        );
        rTree_.remove(std::cref(item));
    }
}

template<class T>
void Group<T>::remove(const Group<T> &other)
{
    auto itr = std::remove_if(items_.begin(), items_.end(), [&other, this](const T &item) -> bool
    {
        if(other.isInGroup(item))
        {
            itemSet_.erase(std::cref(item));
            rTree_.remove(std::cref(item));
            return true;
        }

        return false;
    });

    items_.erase(itr, items_.end());
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
    if (this != &rhs)
        remove(rhs);
    else
        clear();

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
    namespace bgi = boost::geometry::index;

    std::vector<Ref<const T>> result;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
        {
            const Circle &c = static_cast<const Circle &>(shape);
            auto covered = [&c](const T &item) { return c.isInside(item.centroid()); };
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


    return result;
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsCoveredBy(const Shape2D &shape) const
{
    std::vector<Ref<const T>> result;
    result.reserve(items_.size());

    namespace bgi = boost::geometry::index;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
        {
            const Circle &c = static_cast<const Circle &>(shape);
            auto covered = [&c](const T &item) { return c.isCovered(item.centroid()); };
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

    return result;
}

template<class T>
std::vector<Ref<const T> > Group<T>::nearestItems(const Point2D &pt, size_t k) const
{
    std::vector<Ref<const T>> result;
    result.reserve(k);

    rTree_.query(boost::geometry::index::nearest(pt, k),
                 std::back_inserter(result));

    return result;
}

template<class T>
std::vector<Ref<const T>> Group<T>::nearestItems(const Shape2D &shape, size_t k) const
{
    std::vector<Ref<const T>> result;
    result.reserve(k);

    switch (shape.type())
    {
        case Shape2D::CIRCLE:

            rTree_.query(boost::geometry::index::nearest(shape.centroid(), k),
                         std::back_inserter(result));

            break;

        default:
            throw Exception("Group<T>", "nearestItems", "shape type not supported.");
    }

    return result;
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
    return itemSet_.find(std::cref(item)) != itemSet_.end();
}

//- Private helper methods

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