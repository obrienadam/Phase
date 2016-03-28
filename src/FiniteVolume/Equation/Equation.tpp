#include "Equation.h"
#include "Exception.h"

template<class T>
Scalar Equation<T>::solve()
{
    field_ = spMat_.solve(b_);
    return error();
}

template<class T>
Equation<T>& Equation<T>::operator =(const Term& term)
{
    spMat_.assemble(term.coefficients());
    const auto &sources = term.sources();

    for(int i = 0, end = b_.size(); i < end; ++i)
        b_[i] = sources[i];

    return *this;
}
