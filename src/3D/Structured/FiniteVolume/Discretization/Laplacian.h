#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "Structured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Field<T> &phi);

}

#include "Laplacian.tpp"

#endif
