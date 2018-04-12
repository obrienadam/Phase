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
    items_.reserve(size() + items.size());
    itemSet_.reserve(size() + items.size());

    for (const T &item: items)
        if(itemSet_.insert(std::cref(item)).second)
            items_.push_back(std::cref(item));

    rTree_.insert(items.items_.begin(), items.items_.end());
}

template<class T>
void Group<T>::remove(const T &item)
{
    if (itemSet_.erase(std::cref(item)))
    {
        items_.erase(
                std::find_if(items_.begin(), items_.end(), [&item](const T &arg)
                { return item.id() == arg.id(); })
        );

        rTree_.remove(std::cref(item));
    }
}

template<class T>
void Group<T>::remove(const Group<T> &other)
{
    auto itr = std::remove_if(items_.begin(), items_.end(), [&other, this](const T &item) -> bool
    {
        if (other.isInGroup(item))
        {
            itemSet_.erase(std::cref(item));
            return true;
        }

        return false;
    });

    items_.erase(itr, items_.end());
    rTree_.remove(other.begin(), other.end());
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

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
            return std::vector<Ref<const T>>(
                    rTree_.qbegin(bgi::within(shape.boundingBox())
                                  && bgi::satisfies([&shape](const T &item)
                                                    { return shape.isInside(item.centroid()); })),
                    rTree_.qend());
        case Shape2D::BOX:
            return std::vector<Ref<const T>>(rTree_.qbegin(bgi::within(shape.boundingBox())), rTree_.qend());
        case Shape2D::POLYGON:
            return std::vector<Ref<const T>>(rTree_.qbegin(bgi::within(static_cast<const Polygon &>(shape).boostRing())),
                                      rTree_.qend());
    }
}

template<class T>
std::vector<Ref<const T> > Group<T>::itemsCoveredBy(const Shape2D &shape) const
{
    namespace bgi = boost::geometry::index;

    switch (shape.type())
    {
        case Shape2D::CIRCLE:
            return std::vector<Ref<const T>>(
                    rTree_.qbegin(bgi::covered_by(shape.boundingBox())
                                  && bgi::satisfies([&shape](const T &item)
                                                    { return shape.isCovered(item.centroid()); })),
                    rTree_.qend());

        case Shape2D::BOX:
            return std::vector<Ref<const T>>(rTree_.qbegin(bgi::covered_by(shape.boundingBox())), rTree_.qend());
        case Shape2D::POLYGON:
            return std::vector<Ref<const T>>(rTree_.qbegin(bgi::covered_by(static_cast<const Polygon &>(shape).boostRing())),
                                      rTree_.qend());
    }
}

template<class T>
std::vector<Ref<const T> > Group<T>::nearestItems(const Point2D &pt, size_t k) const
{
    namespace bgi = boost::geometry::index;
    return std::vector<Ref<const T>>(rTree_.qbegin(bgi::nearest(pt, k)), rTree_.qend());
}

template<class T>
const T &Group<T>::nearestItem(const Point2D &pt) const
{
    return nearestItems(pt, 1)[0];
}

template<class T>
bool Group<T>::isInGroup(const T &item) const
{
    return itemSet_.find(std::cref(item)) != itemSet_.end();
}

template<class T>
std::vector<Point2D> Group<T>::coordinates() const
{
    std::vector<Point2D> coords(items_.size());
    std::transform(items_.begin(), items_.end(), coords.begin(), [](const T &item)
    { return item.centroid(); });
    return coords;
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