#include "Group.h"

template<class T>
Group<T> &Group<T>::operator=(const std::vector<Ref<const T> > &rhs)
{
    Set<T>::operator =(rhs);
    return *this;
}

template<class T>
void Group<T>::clear()
{
    Set<T>::clear();
    rTree_.clear();
}

template<class T>
bool Group<T>::add(const T &item)
{
    if (Set<T>::add(item))
    {
        rTree_.insert(std::cref(item));
        return true;
    }

    return false;
}

template<class T>
void Group<T>::add(typename Set<T>::const_iterator begin, typename Set<T>::const_iterator end)
{
    Set<T>::add(begin, end);
    rTree_.insert(begin, end);
}

template<class T>
void Group<T>::add(typename std::vector<T>::const_iterator begin, typename std::vector<T>::const_iterator end)
{
    Set<T>::add(begin, end);
    rTree_.insert(begin, end);
}

template<class T>
void Group<T>::add(const Set<T> &set)
{
    Set<T>::add(set);
    rTree_.insert(set.begin(), set.end());
}

template<class T>
bool Group<T>::remove(const T &item)
{
    if (Set<T>::remove(item))
    {
        rTree_.remove(std::cref(item));
        return true;
    }

    return false;
}

template<class T>
void Group<T>::remove(typename Set<T>::const_iterator begin, typename Set<T>::const_iterator end)
{
    Set<T>::remove(begin, end);
    rTree_.remove(begin, end);
}

template<class T>
void Group<T>::remove(const Set<T> &set)
{
    Set<T>::remove(set);
    rTree_.remove(set.begin(), set.end());
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
std::vector<Point2D> Group<T>::coordinates() const
{
    std::vector<Point2D> coords(Set<T>::size());
    std::transform(Set<T>::begin(), Set<T>::end(), coords.begin(), [](const T &item)
    { return item.centroid(); });
    return coords;
}
