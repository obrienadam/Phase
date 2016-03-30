#include <algorithm>

#include "Field.h"

template <class T>
Field<T>::Field(size_t size, const T &initialValue, const std::string &name)
    :
      std::vector<T>(size, initialValue),
      name(name)
{

}

template <class T>
void Field<T>::save()
{
    prev_ = std::vector<T>(*this);
}

template <class T>
bool Field<T>::previousFieldAvailable() const
{
    return prev_.size() == this->size() ? true : false;
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
