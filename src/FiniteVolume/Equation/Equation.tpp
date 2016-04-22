#include <stdio.h>

#include "Equation.h"
#include "Exception.h"

template<class T>
Equation<T>::Equation(const Equation<T>& other)
    :
      spMat_(other.spMat_),
      boundaries_(other.boundaries_),
      sources_(other.sources_),
      field_(other.field_)
{

}

template<class T>
Scalar Equation<T>::solve()
{
    field_ = spMat_.solve(boundaries_ + sources_);
    printf("Solved %s equation, error = %lf, number of iterations = %d\n", name.c_str(), error(), iterations());
    return error();
}

template<class T>
Equation<T>& Equation<T>::operator +=(const Equation<T>& rhs)
{
    spMat_ += rhs.spMat_;
    boundaries_ += rhs.boundaries_;
    sources_ += rhs.sources_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator -=(const Equation<T>& rhs)
{
    spMat_ -= rhs.spMat_;
    boundaries_ -= rhs.boundaries_;
    sources_ -= rhs.sources_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator =(const Equation<T>& rhs)
{
    if(this != &rhs)
    {
        spMat_ = rhs.spMat_;
        boundaries_ = rhs.boundaries_;
        sources_ = rhs.sources_;
    }

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator *=(Scalar rhs)
{
    spMat_ *= rhs;
    boundaries_ *= rhs;
    sources_ *= rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(Scalar rhs)
{
    for(int i = 0, end = boundaries_.rows(); i < end; ++i)
        sources_(i) += rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(const Equation<T>& rhs)
{
    spMat_ -= rhs.spMat_;
    boundaries_ -= rhs.boundaries_;
    sources_ = rhs.sources_ - sources_;
    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator ==(const T& rhs)
{
    for(const Cell& cell: rhs.grid.activeCells())
        sources_(cell.globalIndex()) += rhs[cell.globalIndex()];

    return *this;
}

//- External functions

template<class T>
Equation<T> operator+(Equation<T> lhs, const Equation<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const Equation<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(Equation<T> lhs, const T& rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
Equation<T> operator+(const T& lhs, Equation<T> rhs)
{
    rhs += lhs;
    return rhs;
}

template<class T>
Equation<T> operator-(Equation<T> lhs, const T& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
Equation<T> operator-(const T& lhs, Equation<T> rhs)
{
    rhs -= lhs;
    return rhs;
}

template<class T>
Equation<T> operator*(Equation<T> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class T>
Equation<T> operator*(Scalar lhs, Equation<T> rhs)
{
    rhs *= lhs;
    return rhs;
}
