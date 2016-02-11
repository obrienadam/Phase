#include "Equation.h"

template<class T>
void Equation<T>::solve()
{
    field_ = spMat_.solve(b_);
}
