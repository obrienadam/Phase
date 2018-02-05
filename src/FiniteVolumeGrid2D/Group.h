#ifndef GROUP_H
#define GROUP_H

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "Shape2D.h"
#include "Circle.h"
#include "Box.h"

template<class T>
class Group
{
public:

    typedef typename std::vector<Ref<const T> >::iterator iterator;
    typedef typename std::vector<Ref<const T> >::const_iterator const_iterator;
    typedef boost::geometry::index::quadratic<8, 4> Parameters;

    struct IndexableGetter
    {
        typedef Point2D result_type;

        result_type operator()(const T &item) const
        {
            return item.centroid();
        }
    };

    struct EqualTo
    {
        bool operator()(const T &lhs, const T &rhs) const
        {
            return lhs.id() == rhs.id();
        }
    };

    struct Hash
    {
        std::size_t operator()(const T &item) const
        {
            return std::hash<Label>()(item.id());
        }
    };

    Group(const std::string &name = "N/A") : name_(name)
    {}

    //- Group name
    const std::string &name() const
    { return name_; }

    //- Sizing
    virtual void clear();

    size_t size() const
    { return items_.size(); }

    bool empty() const
    { return items_.empty(); }

    //- Data access
    const std::vector<Ref<const T> > &items() const
    { return items_; }

    //- Adding/removing items
    virtual void add(const T &item);

    virtual void add(const Group<T> &items);

    template<class iterator>
    void add(iterator first, iterator last)
    {
        for (; first != last; ++first)
            add(*first);
    }

    virtual void remove(const T &item);

    virtual void remove(const Group<T> &other);

    template<class iterator>
    void remove(iterator begin, iterator end)
    {
        std::unordered_set<Ref<const T>, Hash, EqualTo> items(begin, end);

        auto itr = std::remove_if(items_.begin(), items_.end(), [&items, this](const T &item)
        {
            if(items.find(std::cref(item)) != items.end())
            {
                itemSet_.erase(std::cref(item));
                rTree_.remove(std::cref(item));
                return true;
            }

            return false;
        });

        items_.erase(itr, items_.end());
    }

    //- Operators
    Group<T> &operator+=(const Group<T> &rhs);

    Group<T> &operator-=(const Group<T> &rhs);

    Group<T> intersection(const Group<T> &other) const;

    Group<T> difference(const Group<T> &other) const;

    Group<T> symmetricDifference(const Group<T> &other) const;

    //- Searching
    std::vector<Ref<const T> > itemsWithin(const Shape2D &shape) const;

    std::vector<Ref<const T> > itemsNotWithin(const Shape2D &shape) const;

    std::vector<Ref<const T> > itemsCoveredBy(const Shape2D &shape) const;

    std::vector<Ref<const T> > itemsNotCoveredBy(const Shape2D &shape) const;

    std::vector<Ref<const T> > nearestItems(const Point2D &pt, size_t k) const;

    std::vector<Ref<const T> > nearestItems(const Shape2D &shape, size_t k) const;

    const T &nearestItem(const Point2D &pt) const;

    const T &nearestItem(const Shape2D &shape) const;

    //- Iterators
    iterator begin()
    { return items_.begin(); }

    iterator end()
    { return items_.end(); }

    const_iterator begin() const
    { return items_.begin(); }

    const_iterator end() const
    { return items_.end(); }

    //- Access
    const std::unordered_set<Ref<const T>, Hash, EqualTo> &itemSet() const
    { return itemSet_; }

    bool isInGroup(const T &item) const;

protected:

    //- Data
    std::string name_;
    std::unordered_set<Ref<const T>, Hash, EqualTo> itemSet_; // Allows cell lookup via an id
    std::vector<Ref<const T> > items_; // Used for faster iteration over all cells
    boost::geometry::index::rtree<Ref<const T>, Parameters, IndexableGetter, EqualTo> rTree_; //- For searching
};

template<class T>
Group<T> operator+(Group<T> lhs, const Group<T> &rhs);

template<class T>
Group<T> operator-(Group<T> lhs, const Group<T> &rhs);

#include "Group.tpp"

#endif
