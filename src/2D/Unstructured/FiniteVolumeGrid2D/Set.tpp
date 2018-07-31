#include "Set.h"

template<class T>
void Set<T>::clear()
{
    items_.clear();
    itemSet_.clear();
}

template<class T>
Set<T>& Set<T>::operator+=(const Set &rhs)
{
    if(this != &rhs)
        add(rhs.begin(), rhs.end());
    return *this;
}

template<class T>
Set<T> &Set<T>::operator-=(const Set &rhs)
{
    remove(rhs);
    return *this;
}

template<class T>
Set<T> &Set<T>::operator=(const std::vector<Ref<const T>> &rhs)
{
    clear();
    add(rhs.begin(), rhs.end());
    return *this;
}

template<class T>
void Set<T>::reserve(Size size)
{
    items_.reserve(size);
    itemSet_.reserve(size);
}

template<class T>
bool Set<T>::add(const T &item)
{
    auto insert = itemSet_.insert(std::cref(item));

    if(insert.second)
        items_.push_back(std::cref(item));

    return insert.second;
}

template<class T>
void Set<T>::add(const_iterator begin, const_iterator end)
{
    items_.reserve(size() + (end - begin));
    itemSet_.reserve(size() + (end - begin));

    for(auto itr = begin; itr != end; ++itr)
        if(itemSet_.insert(*itr).second)
            items_.push_back(*itr);
}

template<class T>
void Set<T>::add(typename std::vector<T>::const_iterator begin, typename std::vector<T>::const_iterator end)
{
    items_.reserve(size() + (end - begin));
    itemSet_.reserve(size() + (end - begin));

    for(auto itr = begin; itr != end; ++itr)
        if(itemSet_.insert(*itr).second)
            items_.push_back(*itr);
}

template<class T>
void Set<T>::add(const Set<T> &set)
{
    add(set.begin(), set.end());
}

template<class T>
bool Set<T>::remove(const T &item)
{
    bool removed = itemSet_.erase(std::cref(item));

    if(removed)
    {
        items_.erase(
                    std::find_if(items_.begin(), items_.end(), [&item](const T &arg)
        { return item.id() == arg.id(); })
                    );
    }

    return removed;
}

template<class T>
void Set<T>::remove(const_iterator begin, const_iterator end)
{
    std::unordered_set<Ref<const T>, Hash, EqualTo> items(begin, end);

    auto itr = std::remove_if(items_.begin(), items_.end(), [&items, this](const T &item) -> bool
    {
        if (items.find(std::cref(item)) != items.end())
        {
            itemSet_.erase(std::cref(item));
            return true;
        }

        return false;
    });

    items_.erase(itr, items_.end());
}

template<class T>
void Set<T>::remove(const Set<T> &set)
{
    auto itr = std::remove_if(items_.begin(), items_.end(), [&set, this](const T &item) -> bool
    {
        if (set.isInSet(item))
        {
            itemSet_.erase(std::cref(item));
            return true;
        }

        return false;
    });

    items_.erase(itr, items_.end());
}
