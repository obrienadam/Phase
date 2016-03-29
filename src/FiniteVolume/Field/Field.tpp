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
