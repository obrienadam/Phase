#include <algorithm>

#include "Set.h"

template<class T>
void Set<T>::add(const T &item)
{
    if(_itemSet.insert(item).second)
        _items.push_back(item);
}

template<class T>
template<class InputIterator>
void Set<T>::add(InputIterator first, InputIterator last)
{
    _itemSet.reserve(last - first);
    _items.reserve(last - first);

    for(auto it = first; it != last; ++it)
        add(*it);
}

template<class T>
void Set<T>::remove(const T &item)
{
    if((bool)_itemSet.erase(item))
        _items.erase(std::find_if(_items.begin(), _items.end(), [&item](const T &val)
        { return item.id() == val.id(); }));
}

template<class T>
template<class InputIterator>
void Set<T>::remove(InputIterator first, InputIterator last)
{
    for(auto it = first; it != last; ++it)
        _itemSet.erase(*it);

    _items.erase(std::remove_if(_items.begin(), _items.end(), [this](const T &item)
    {
        return _itemSet.find(item) == _itemSet.end();
    }), _items.end());
}
