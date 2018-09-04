#ifndef TYPES_H
#define TYPES_H
// This file can be used to alter the working precision of the code

#include <functional>

typedef double Scalar;
typedef std::size_t Label;
typedef std::size_t Size;
typedef int Index;

template <class T>
using Ref = std::reference_wrapper<T>;

#endif
