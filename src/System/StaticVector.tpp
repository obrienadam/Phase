#include "StaticVector.h"

template<class T, std::size_t N>
typename StaticVector<T, N>::iterator StaticVector<T, N>::insert(iterator pos, const T &value)
{
    ++_size;
    std::copy(rbegin() + 1, rbegin() + _size - (pos - begin()), rbegin());
    *pos = value;
    return pos;
}

template<class T, std::size_t N>
template<class InputIt>
typename StaticVector<T, N>::iterator StaticVector<T, N>::insert(iterator pos,
                                                                 InputIt first,
                                                                 InputIt last)
{
    _size += last - first;
    std::copy(rbegin() + (last - first), rbegin() + _size - (pos - begin()), rbegin());
    std::copy(first, last, pos);
    return pos;
}
