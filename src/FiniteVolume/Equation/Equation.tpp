#include "Equation.h"
#include "Exception.h"

template<class T>
Scalar Equation<T>::solve()
{
    field_ = spMat_.solve(b_);
    return error();
}
