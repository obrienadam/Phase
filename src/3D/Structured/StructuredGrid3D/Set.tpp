#include <algorithm>

#include "Set.h"

template<class T>
bool Set<T>::add(const T& item)
{
    auto insert = _itemSet.insert(std::cref(item));
    if(insert.second)
        _items.push_back(std::cref(item));
    return insert.second;
}

template<class T>
bool Set<T>::remove(const T& item)
{
    bool removed = _itemSet.erase(std::cref(item));
    if(removed)
        _items.erase(std::find_if(_items.begin(), _items.end(), [&item](const T &cmp) {
            return cmp.id() == item.id();
        }));

    return removed;
}
