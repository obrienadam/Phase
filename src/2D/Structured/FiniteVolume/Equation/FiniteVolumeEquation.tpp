#include "FiniteVolumeEquation.h"

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(FiniteVolumeEquation<T> &&rhs)
{
    Equation::operator =(rhs);
    return *this;
}
