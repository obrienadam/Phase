#include <algorithm>

#include "Field.h"

template <class T>
Field<T>::Field(int sizeI, int sizeJ, const T &initialValue, const std::string &name)
    :
      sizeI_(sizeI),
      sizeJ_(sizeJ),
      std::vector<T>(sizeI*sizeJ, initialValue),
      name(name)
{

}

template <class T>
Field<T>::Field(int sizeI, int sizeJ, const std::string &name)
    :
      Field(sizeI, sizeJ, T(), name)
{

}

template <class T>
T& Field<T>::operator ()(int i, int j)
{
    return std::vector<T>::operator [](i + sizeI_*j);
}

template <class T>
const T& Field<T>::operator ()(int i, int j) const
{
    return std::vector<T>::operator [](i + sizeI_*j);
}

template <class T>
void Field<T>::fill(const T &value)
{
    std::fill(this->begin(), this->end(), value);
}

template <class T>
void Field<T>::resize(int sizeI, int sizeJ)
{
    sizeI_ = sizeI;
    sizeJ_ = sizeJ;
    std::vector<T>::resize(sizeI_*sizeJ_);
}
