#ifndef PHASE_STATIC_VECTOR_H
#define PHASE_STATIC_VECTOR_H

#include <array>

template<class T, std::size_t N>
class StaticVector
{
public:

    typedef typename std::array<T, N>::iterator iterator;

    typedef typename std::array<T, N>::const_iterator const_iterator;

    StaticVector(): _begin(_data.begin()), _end(_data.begin()) {}

    void push_back(const T& val)
    { *(_end++) = val; }

    std::size_t size() const
    { return _end - _begin; }

    std::size_t capacity() const
    { return _data.size(); }

    //- Iterators

    iterator begin()
    { return _begin; }

    iterator end()
    { return _end; }

    const_iterator begin() const
    { return _begin; }

    const_iterator end() const
    { return _end; }

private:

    std::array<T, N> _data;

    typename std::array<T, N>::iterator _begin, _end;
};

#include "StaticVector.tpp"

#endif
