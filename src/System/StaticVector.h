#ifndef PHASE_STATIC_VECTOR_H
#define PHASE_STATIC_VECTOR_H

#include <array>

template<class T, std::size_t N>
class StaticVector
{
public:

    typedef typename std::array<T, N>::iterator iterator;

    typedef typename std::array<T, N>::const_iterator const_iterator;

    typedef typename std::array<T, N>::reverse_iterator reverse_iterator;

    typedef typename std::array<T, N>::const_reverse_iterator const_reverse_iterator;

    StaticVector(std::size_t size = 0) : _size(size) {}

    StaticVector(std::size_t size, const T &val) : _size(size)
    { std::fill_n(_data.begin(), _size, val); }

    //- Operators
    T& operator[](std::size_t i)
    { return _data[i]; }

    const T& operator [](std::size_t i) const
    { return _data[i]; }

    void push_back(const T& val)
    { *(_data.begin() + _size++) = val; }

    iterator insert(iterator pos, const T &value);

    template<class InputIt>
    iterator insert(iterator pos, InputIt first, InputIt last);

    iterator insert(iterator pos, const std::initializer_list<T> &ilist)
    { return insert(pos, ilist.begin(), ilist.end()); }

    std::size_t size() const
    { return _size; }

    std::size_t capacity() const
    { return N; }

    bool empty() const
    { return (bool)_size; }

    void resize(std::size_t size)
    { _size = size; }

    void pop_back()
    { _size = _size > 0 ? _size - 1: 0; }

    void clear()
    { _size = 0; }

    //- Iterators

    iterator begin()
    { return _data.begin(); }

    iterator end()
    { return _data.begin() + _size; }

    const_iterator begin() const
    { return _data.begin(); }

    const_iterator end() const
    { return _data.begin() + _size; }

    reverse_iterator rbegin()
    { return _data.rbegin() + N - _size; }

    reverse_iterator rend()
    { return _data.rend(); }

    const_iterator cbegin() const
    { return _data.cbegin(); }

    const_iterator cend() const
    { return _data.cbegin() + _size; }

private:

    std::array<T, N> _data;

    std::size_t _size;
};

#include "StaticVector.tpp"

#endif
