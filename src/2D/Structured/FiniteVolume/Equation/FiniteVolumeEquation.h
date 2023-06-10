#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include "FiniteVolume/Field/Field.h"

#include "Math/Equation.h"

template <class T> class FiniteVolumeEquation : public Equation {
public:
  //- Constructors
  FiniteVolumeEquation(Field<T> &field, int nnz = 5);

  FiniteVolumeEquation(const std::string &name, Field<T> &field, int nnz = 5)
      : FiniteVolumeEquation(field, nnz) {
    _name = name;
  }

  FiniteVolumeEquation(const FiniteVolumeEquation<T> &other) = default;

  FiniteVolumeEquation(FiniteVolumeEquation<T> &&other) = default;

  //- Operators
  FiniteVolumeEquation &operator=(FiniteVolumeEquation<T> &&rhs);

  //- Add/remove coefficients

  void add(const Cell &cell0, const Cell &cell1, Scalar coeff);

  void addSource(const Cell &cell, const T &src);

  //- name access
  const std::string &name() const { return _name; }

private:
  std::string _name;

  Field<T> &_field;
};

#include "FiniteVolumeEquation.tpp"

#endif
