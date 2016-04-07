#include "Equation.h"
#include "Exception.h"

template<class T>
Equation<T>::Equation(const Equation<T>& other)
    :
      spMat_(other.spMat_),
      b_(other.b_),
      field_(other.field_)
{

}

template<class T>
Scalar Equation<T>::solve()
{
    field_ = spMat_.solve(b_);
    return error();
}

template<class T>
Equation<T>& Equation<T>::operator +=(const Equation<T>& rhs)
{
    spMat_ += rhs.spMat_;
    b_ += rhs.b_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator -=(const Equation<T>& rhs)
{
    spMat_ -= rhs.spMat_;
    b_ -= rhs.b_;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator =(const Equation<T>& rhs)
{
    if(this != &rhs)
    {
        spMat_ = rhs.spMat_;
        b_ = rhs.b_;
    }

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator *=(Scalar rhs)
{
    spMat_ *= rhs;
    b_ *= rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(Scalar rhs)
{
    for(int i = 0, end = b_.rows(); i < end; ++i)
        b_[i] += rhs;

    return *this;
}

template<class T>
Equation<T>& Equation<T>::operator==(const Equation<T>& rhs)
{
    return operator+=(rhs);
}

template<class T>
Equation<T>& Equation<T>::operator ==(const T& rhs)
{
    return Equation<T>::operator -=(rhs);
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
