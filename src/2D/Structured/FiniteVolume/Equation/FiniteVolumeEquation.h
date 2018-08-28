#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include "Math/Equation.h"
#include "FiniteVolume/Field/Field.h"

template<class T>
class FiniteVolumeEquation: public Equation
{
public:

    FiniteVolumeEquation(Field<T> &field);

private:

    Field<T> &_field;
};

#endif
