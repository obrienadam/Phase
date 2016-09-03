#include <algorithm>

#include "Field.h"

template <class T>
Field<T>::Field(size_t size, const T &initialValue, const std::string &name)
    :
      std::vector<T>(size, initialValue),
      name_(name)
{

}

//- External functions
