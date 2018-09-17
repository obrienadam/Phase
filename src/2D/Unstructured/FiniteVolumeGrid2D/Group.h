#ifndef PHASE_GROUP_H
#define PHASE_GROUP_H

#include <string>
#include <unordered_set>

#include "Geometry/Shape2D.h"
#include "Geometry/Circle.h"
#include "Geometry/Box.h"

#include "Set.h"

template<class T>
class Group: public Set<T>
{
public:

    typedef boost::geometry::index::quadratic<8, 4> Parameters;

    struct IndexableGetter
    {
        typedef Point2D result_type;

        result_type operator()(const T &item) const
        { return item.centroid(); }
    };

    Group(const std::string &name = "") : Set<T>(name)
    {}

    template<class const_iterator>
    Group(const_iterator first, const_iterator last, const std::string &name = "") : Group(name)
    { add(first, last); }

    //- Assignment
    Group<T> &operator=(const std::vector<Ref<const T>> &rhs);

    //- Sizing
    void clear() override;

    //- Adding/removing items
    virtual bool add(const T &item) override;

    virtual void add(typename Set<T>::const_iterator begin, typename Set<T>::const_iterator end) override;

    virtual void add(typename std::vector<T>::const_iterator begin, typename std::vector<T>::const_iterator end) override;

    virtual void add(const Set<T> &set) override;

    virtual bool remove(const T &item) override;

    virtual void remove(typename Set<T>::const_iterator begin, typename Set<T>::const_iterator end) override;

    virtual void remove(const Set<T> &set) override;

    //- Searching
    std::vector<Ref<const T> > itemsWithin(const Shape2D &shape) const;

    std::vector<Ref<const T> > itemsCoveredBy(const Shape2D &shape) const;

    void itemsWithin(const Shape2D &shape, std::vector<Ref<const T>> &result) const;

    void itemsCoveredBy(const Shape2D &shape, std::vector<Ref<const T>> &result) const;

    std::vector<Ref<const T> > nearestItems(const Point2D &pt, size_t k) const;

    const T &nearestItem(const Point2D &pt) const;

    //- Misc

    std::vector<Point2D> coordinates() const;

protected:

    boost::geometry::index::rtree<Ref<const T>, Parameters, IndexableGetter, typename Set<T>::EqualTo> rTree_; //- For searching
};

#include "Group.tpp"

#endif
