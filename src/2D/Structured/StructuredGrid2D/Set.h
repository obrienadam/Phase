#ifndef PHASE_SET_H
#define PHASE_SET_H

#include <unordered_set>
#include <vector>

#include "Types/Types.h"

template<class T>
class Set
{
public:

    typedef typename std::vector<Ref<const T>>::iterator iterator;

    typedef typename std::vector<Ref<const T>>::const_iterator const_iterator;

    Set(const std::string &name = "") : _name(name)
    { }

    //- insertion/deletion

    bool add(const T &item);

    template<class InputIterator>
    void add(InputIterator first, InputIterator last);

    void remove(const T &item);

    template<class InputIterator>
    void remove(InputIterator first, InputIterator last);

    std::size_t size() const
    { return _items.size(); }

    //- Iterators

    iterator begin()
    { return _items.begin(); }

    iterator end()
    { return _items.end(); }

    const_iterator begin() const
    { return _items.begin(); }

    const_iterator end() const
    { return _items.end(); }

    //- Tests
    bool isInSet(const T& item) const
    { return _itemSet.find(item) != _itemSet.end(); }

protected:

    struct Hash
    {
        std::size_t operator()(const T &item) const
        { return std::hash<std::size_t>{}(item.id()); }
    };

    struct EqualTo
    {
        bool operator()(const T &lhs, const T &rhs) const
        { return lhs.id() == rhs.id(); }
    };

    std::string _name;

    std::unordered_set<Ref<const T>, Hash, EqualTo> _itemSet;

    std::vector<Ref<const T>> _items;
};

#include "Set.tpp"

#endif
