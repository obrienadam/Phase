#ifndef GROUP_H
#define GROUP_H

#include <string>
#include <unordered_map>

#include "Shape2D.h"
#include "Circle.h"
#include "Box.h"

template <class T>
class Group
{
public:

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
    std::vector<Ref<const T> > items() const
    { return items_; }

    //- Adding/removing items
    virtual void add(const T &item);

    virtual void add(const Group<T>& other);

    template<class const_iterator>
    void add(const_iterator begin, const_iterator end)
    {
        for(; begin != end; ++begin)
            add(*begin);
    }

    virtual void remove(const T &item);

    virtual void remove(const Group<T>& other);

    virtual void merge(Group<T> &other);

    //- Searching
    std::vector<Ref<const T> > itemsWithin(const Shape2D &shape) const;

    std::vector<Ref<const T> > itemsWithin(const Circle &circle) const;

    std::vector<Ref<const T> > itemsWithin(const Box &box) const;

    std::vector<Ref<const T> > nearestItems(const Point2D &pt, size_t k) const;

    //- Iterators
    typename std::vector<Ref<const T> >::iterator begin()
    { return items_.begin(); }

    typename std::vector<Ref<const T> >::iterator end()
    { return items_.end(); }

    typename std::vector<Ref<const T> >::const_iterator begin() const
    { return items_.begin(); }

    typename std::vector<Ref<const T> >::const_iterator end() const
    { return items_.end(); }

    bool isInGroup(const T &item) const;

protected:

    typedef std::pair<Point2D, Label> Value;
    typedef boost::geometry::index::rtree<Value, boost::geometry::index::quadratic<16> > Rtree;

    std::vector<Ref<const T> > getRefs(const std::vector<Value> &vals) const;

    std::string name_;
    std::unordered_map<Label, Ref<const T> > itemSet_; // Allows cell lookup via an id
    std::vector<Ref<const T> > items_; // Used for faster iteration over all cells

    Rtree rTree_; //- For searching
};

#include "Group.tpp"

#endif
