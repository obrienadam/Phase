#include <algorithm>

#include "Field.h"

template <class T>
Field<T>::Field(size_t size, const T &initialValue, const std::string &name)
    :
      std::vector<T>(size, initialValue),
      name(name),
      prevFieldPtr_(nullptr)
{

}

template <class T>
void Field<T>::save()
{
    prevFieldPtr_ = std::shared_ptr< Field<T> >(new Field<T>(*this));
}

template <class T>
Field<T>& Field<T>::operator +=(const Field<T>& rhs)
{
    auto& self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] += rhs[i];

    return self;
}

template <class T>
Field<T>& Field<T>::operator -=(const Field<T>& rhs)
{
    auto& self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] -= rhs[i];

    return self;
}

template <class T>
template <class T2>
Field<T>& Field<T>::operator *=(const Field<T2>& rhs)
{
    auto& self = *this;

    for(int i = 0, end = self.size(); i < end; ++i)
        self[i] *= rhs[i];

    return self;
}

template<class T>
T Field<T>::max() const
{
    T maxVal();

    for(const T& val: *this)
        maxVal = val > maxVal ? val : maxVal;

    return maxVal;
}

//- External functions

template <class T>
Field<T> operator+(Field<T> lhs, const Field<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
Field<T> operator-(Field<T> lhs, const Field<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}
