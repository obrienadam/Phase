#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include "Math/Equation.h"
#include "FiniteVolume/Field/Field.h"

template<class T>
class FiniteVolumeEquation: public Equation
{
public:

    //- Constructors
    FiniteVolumeEquation(Field<T> &field);

    FiniteVolumeEquation(const std::string &name, Field<T> &field) : FiniteVolumeEquation(field)
    { _name = name; }

    FiniteVolumeEquation(const FiniteVolumeEquation<T> &other) = default;

    FiniteVolumeEquation(FiniteVolumeEquation<T> &&other) = default;

    //- Operators
    FiniteVolumeEquation &operator =(FiniteVolumeEquation<T> &&rhs);

    //- name access
    const std::string &name() const
    { return _name; }

private:

    std::string _name;

    Field<T> &_field;
};

#include "FiniteVolumeEquation.tpp"

#endif
